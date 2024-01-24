# Methods to compute the multi-order Laplacian of a hypergraph.
# Directed hypergraphs are treated as undirected.
import sys
sys.path.insert(1, '../')
from helpers.io import read_as_undirected_hyperedge_list
from helpers.parser_for_output_file_names import parser_for_output_file_name as parse_file_name
import zipfile
import io
from datetime import datetime
from scipy.sparse import lil_matrix, csr_matrix, coo_matrix
import numpy as np
from sklearn.utils import check_array, check_symmetric
from scipy import sparse
from scipy.sparse.csgraph import laplacian as csgraph_laplacian
from scipy.sparse.linalg import lobpcg
from pyamg import smoothed_aggregation_solver
from itertools import combinations


def get_adj_matrix(edges, N):
    # creates the adjacency matrix of the hypergraph
    # given its hyperedges
    A = lil_matrix((N, N))
    for e in edges:
        for v1, v2 in combinations(e, 2):
            A[v1, v2] += 1
            A[v2, v1] += 1
    return csr_matrix(A)


def compute_eigenvalues(edges, max_d, eig_num, tol, N, gamma, symmetric):
    '''
    INPUT
    ======
    edges, list: undirected hyperedges
    max_d, int: max order of laplacian to compute
    eig_num, int: number of eigenvalues to compute
    tol, float: tolerance when finding the eigenvalues
    N, int: number of vertices in the hypergraph
    gamma, np.array (max_d,): importance of each order
    symmetric, bool: wheteher the Laplacian should be normalized

    OUTPUT
    ======
    eivals, list: eigenvalues of the laplacian of each order d
    eivals_mul, np.array (N,): eigenvalues of the multi-order laplacian
    '''
    eivals = []
    eivals_mul = np.zeros(eig_num)
    L_mul = csr_matrix((N, N))
    print('max_d', max_d)
    for d in range(2, max_d + 1):
        print('examining order d=', d)
        # all edges of size == d
        sub_edges = list(set([e for e in edges if len(e) == d]))
        A = get_adj_matrix(sub_edges, N)
        if len(sub_edges) == 0 or A.count_nonzero() == 0:
            eivals.append([])
            continue
        print("Non-zero", A.count_nonzero(), "Dim", A.shape)
        # d-order eigenvalues
        print(f'computing eigenvalues {datetime.now().strftime("%H:%M:%S")}')
        # we do not multiply the diagonal by d - 1 at this point, but we
        # add d - 1 afterwards, because the eigenvalues shift by c if the
        # diagonal is multiplied by d
        Ld, Kd, eival = find_eigenvalues(adjacency=A,
                                         n_components=min(eig_num, N),
                                         norm_laplacian=symmetric,
                                         eigen_tol=tol,
                                         is_laplacian=False)
        # eival += (d - 1)
        print(f'eigenvalues computed {datetime.now().strftime("%H:%M:%S")}')
        print(eival)
        eivals.append(eival)
        # multi-order laplacian
        if np.sum(Kd) > 0:
            tmp_lap = lil_matrix(Ld)
            for i in range(tmp_lap.shape[0]):
                tmp_lap[i, i] *= (d - 1)
            L_mul += (gamma[d - 2] / np.mean(Kd)) * tmp_lap
    print(
        f'computing multi-order eigenvalues {datetime.now().strftime("%H:%M:%S")}')
    _, _, eivals_mul = find_eigenvalues(adjacency=L_mul,
                                        n_components=min(eig_num, N),
                                        norm_laplacian=False,
                                        eigen_tol=tol,
                                        is_laplacian=True)
    print(
        f'multi-order eigenvalues computed {datetime.now().strftime("%H:%M:%S")}')
    return eivals, eivals_mul


def parallel_laplacian(inp):
    file_path = inp[0]
    file_name = inp[1]
    vmap = inp[2]
    max_d = inp[3]
    N = inp[4]
    gamma = inp[5]
    tol = inp[6]
    eig_num = inp[7]
    symmetric = inp[8]
    typ = inp[9]

    map_fields = parse_file_name(file_name, typ)
    z = zipfile.ZipFile(file_path, "r")
    with z.open(file_name, 'r') as f1:
        fobj = io.TextIOWrapper(f1, encoding='utf-8', newline='')
        uedges = read_as_undirected_hyperedge_list(fobj, vmap)
    eis, eis_m = compute_eigenvalues(uedges,
                                     max_d,
                                     eig_num,
                                     tol,
                                     N,
                                     gamma,
                                     symmetric)

    return [map_fields['algorithm'], map_fields['randomSeed'], eis, eis_m]


def find_eigenvalues(
    adjacency,
    n_components=8,
    seed=0,
    eigen_tol=None,
    norm_laplacian=False,
    is_laplacian=False
):
    """
    Parameters
    ----------
    adjacency : {array-like, sparse graph} of shape (n_samples, n_samples)
        The adjacency matrix of the graph to embed.

    n_components : int, default=8
        The dimension of the projection subspace.

    seed : int,
        A pseudo random number generator used for the initialization
        of the lobpcg eigen vectors decomposition

    eigen_tol : float, default=None
        Stopping criterion for eigendecomposition of the Laplacian matrix.
        Note that values of `tol<1e-5` may lead
        to convergence issues and should be avoided.

    norm_laplacian : bool, default=True
        If True, then compute symmetric normalized Laplacian.
    """

    random_state = np.random.RandomState(seed)
    np.random.seed(seed)

    if not is_laplacian:
        adjacency = check_symmetric(adjacency)
        laplacian, dd = csgraph_laplacian(adjacency,
                                          normed=norm_laplacian,
                                          return_diag=True)
    else:
        laplacian = coo_matrix(adjacency)
        dd = adjacency.diagonal()

    # Use AMG to get a preconditioner and speed up the eigenvalue
    # problem.
    laplacian = check_array(
        laplacian, dtype=[np.float64, np.float32], accept_sparse=True
    )
    laplacian = _set_diag(laplacian, 1, norm_laplacian)

    # The Laplacian matrix is always singular, having at least one zero
    # eigenvalue, corresponding to the trivial eigenvector, which is a
    # constant. Using a singular matrix for preconditioning may result in
    # random failures in LOBPCG and is not supported by the existing
    # theory:
    #     see https://doi.org/10.1007/s10208-015-9297-1
    # Shift the Laplacian so its diagononal is not all ones. The shift
    # does change the eigenpairs however, so we'll feed the shifted
    # matrix to the solver and afterward set it back to the original.
    diag_shift = 1e-5 * sparse.eye(laplacian.shape[0])
    laplacian += diag_shift
    ml = smoothed_aggregation_solver(check_array(laplacian,
                                                 accept_sparse="csr"))
    laplacian -= diag_shift

    M = ml.aspreconditioner()
    # Create initial approximation X to eigenvectors
    X = random_state.standard_normal(size=(laplacian.shape[0],
                                           n_components + 1))
    X[:, 0] = dd.ravel()
    X = X.astype(laplacian.dtype)

    tol = eigen_tol
    eigenvalues, _ = lobpcg(laplacian, X, M=M, tol=tol, largest=False)
    return laplacian, dd, eigenvalues[:n_components]


def _set_diag(laplacian, value, norm_laplacian):
    """Set the diagonal of the laplacian matrix and convert it to a
    sparse format well suited for eigenvalue decomposition.

    Parameters
    ----------
    laplacian : {ndarray, sparse matrix}
        The graph laplacian.

    value : float
        The value of the diagonal.

    norm_laplacian : bool
        Whether the value of the diagonal should be changed or not.

    Returns
    -------
    laplacian : {array, sparse matrix}
        An array of matrix in a form that is well suited to fast
        eigenvalue decomposition, depending on the band width of the
        matrix.
    """
    n_nodes = laplacian.shape[0]
    # We need all entries in the diagonal to values
    if not sparse.isspmatrix(laplacian):
        if norm_laplacian:
            laplacian.flat[:: n_nodes + 1] = value
    else:
        laplacian = laplacian.tocoo()
        if norm_laplacian:
            diag_idx = laplacian.row == laplacian.col
            laplacian.data[diag_idx] = value
        # If the matrix has a small number of diagonals (as in the
        # case of structured matrices coming from images), the
        # dia format might be best suited for matvec products:
        n_diags = np.unique(laplacian.row - laplacian.col).size
        if n_diags <= 7:
            # 3 or less outer diagonals on each side
            laplacian = laplacian.todia()
        else:
            # csr has the fastest matvec and is thus best suited to
            # arpack
            laplacian = laplacian.tocsr()
    return laplacian
