# Compute GENEPY and ECI for directed bipartite graphs representing 
# directed hypergraphs
import sys
sys.path.insert(1, '../')
from helpers.config import folder, side, group, met_path, sample_path, max_workers
from helpers.io import read_as_bipartite
from helpers.parser_for_output_file_names import parser_for_output_file_name as parse_file_name
import numpy as np
import econci
import pandas as pd
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import eigsh
from copy import deepcopy
from tqdm.contrib.concurrent import process_map
import zipfile
import io


def genepy_index(M, group):
    # M, csr_matrix: is the incidence matrix of the bipartite network
    # group, str: indicates which index to
    #             compute between 'countries' and 'products'
    # Node degree k_c
    kc = np.array(M.sum(axis=1)).flatten()
    div = deepcopy(kc)
    div[div == 0] = 1
    # Random walk matrix
    RW = M / div[:, np.newaxis]
    # Compute k_p'
    kp_1 = np.array(RW.sum(axis=0)).flatten()
    # Compute denominator of the fraction
    den = np.outer(kc, kp_1)
    den[den == 0] = 1
    # Compute transformation matrix Wcp
    W = M / den
    # Compute proximity matrix P
    if group == 'countries':
        P = W.dot(W.T)  # P equals Ncc = Wcp*Wcp'
    elif group == 'products':
        P = W.T.dot(W)  # P equals Gpp = Wcp'*Wcp

    # Set diagonal to zero
    P = lil_matrix(P)
    P.setdiag(0)
    P = csr_matrix(P)

    # Compute eigenvectors and eigenvalues of the matrix P
    num_eigenvectors = 2  # Number of eigenvectors to compute
    evalues, evectors = eigsh(P, k=num_eigenvectors, which='LM')

    # Sort absolute eigenvalues in descending order
    order_eig = np.argsort(np.abs(evalues))[::-1]
    # Re-order eigenvalues and eigenvectors accordingly
    evalues = evalues[order_eig]
    evectors = evectors[:, order_eig]
    # First eigenvector (X1 for countries, Y1 for products)
    E1 = np.abs(evectors[:, 0])
    # Second eigenvector (X2 for countries, Y2 for products)
    E2 = evectors[:, 1]

    # Compute GENEPY index using unique contribution
    # Matrix of first two eigenvectors E1 and E2
    E = np.column_stack((E1, E2))
    # Vector of first two eigenvalues
    lambda_ = evalues[:2]
    pow_sum = ((E ** 2) * lambda_).sum(axis=1) ** 2
    GENEPY = pow_sum + 2 * (E ** 2).dot(lambda_ ** 2)
    return E1, E2, GENEPY


def biadjacency_mat(fobj, vmap, side):
    edges = read_as_bipartite(fobj, vmap, side)
    N = len(vmap)
    M = max(x[1] for x in edges) + 1
    print(N, M)

    A = lil_matrix((N, M))
    for e in edges:
        A[e[0], e[1]] = 1
    return csr_matrix(A)


def create_country_product_df(fobj, vmap, side):
    edges = read_as_bipartite(fobj, vmap, side)
    df = pd.DataFrame(edges, columns=['country', 'product'])
    df['export'] = 1
    return df


def parallel_eci(inp):
    '''
    Run ECI for the sample read from file_name
    '''
    file_path = inp[0]
    file_name = inp[1]
    vmap = inp[2]
    side = inp[3]
    typ = inp[4]

    z = zipfile.ZipFile(file_path, "r")
    with z.open(file_name, 'r') as f1:
        map_fields = parse_file_name(file_name, typ)
        fobj = io.TextIOWrapper(f1, encoding='utf-8', newline='')
        df = create_country_product_df(fobj, vmap, side)
    comp = econci.Complexity(df, c='country', p='product', values='export')
    comp.calculate_indexes()
    return comp.eci, map_fields


def parallel_genepy(inp):
    '''
    Run GENEPY for the sample read from file_name
    '''
    file_path = inp[0]
    file_name = inp[1]
    vmap = inp[2]
    side = inp[3]
    group = inp[4]
    typ = inp[5]

    z = zipfile.ZipFile(file_path, "r")
    with z.open(file_name, 'r') as f1:
        map_fields = parse_file_name(file_name, typ)
        fobj = io.TextIOWrapper(f1, encoding='utf-8', newline='')
        M = biadjacency_mat(fobj, vmap, side)
    _, _, genepy = genepy_index(M, group)
    return genepy, map_fields


if __name__ == '__main__':

    root = sys.argv[1]
    typ = sys.argv[2]
    data_path = folder + root
    print(root, typ)

    # run GENEPY observed hypergraph
    vmap = dict()
    with open(f'{data_path}.tsv') as in_f:
        M = biadjacency_mat(in_f, vmap, side)
    _, _, genepy = genepy_index(M, group)
    with open(met_path + f"genepy_{root}_{side}.tsv", 'w') as out_f:
        for idx, v in enumerate(genepy):
            out_f.write(f'{idx}\t{v}\n')
    # run ECI observed hypergraph
    with open(f'{data_path}.tsv') as in_f:
        df = create_country_product_df(in_f, vmap, side)
    comp = econci.Complexity(df, c='country', p='product', values='export')
    comp.calculate_indexes()
    out_file = met_path + f"eci_{root}_{side}.tsv"
    comp.eci.to_csv(out_file, sep="\t", header=False, index=False)
    # save index
    with open(met_path + f"vmap_{root}.tsv", 'w') as o_f:
        for key in vmap:
            o_f.write(f'{key}\t{vmap[key]}\n')

    file_path = f'{sample_path}/{typ}/{root}.zip'
    out_f = open(met_path + f"genepy_{root}_samples_{typ}_{side}.tsv", 'w')
    out_eci = met_path + f"eci_{root}_samples_{typ}_{side}.tsv"

    z = zipfile.ZipFile(file_path, "r")
    zinfo = z.namelist()
    inputs = []
    inputs2 = []
    for file_name in zinfo:
        if file_name.startswith(".") or file_name.startswith("__MACOSX"):
            continue
        if not file_name.endswith('.tsv'):
            continue
        inputs.append([file_path, file_name, vmap, side, group, typ])
        inputs2.append([file_path, file_name, vmap, side, typ])
    # GENEPY
    outputs = process_map(parallel_genepy,
                          inputs,
                          max_workers=max_workers)
    # ECI
    outputs2 = process_map(parallel_eci,
                           inputs2,
                           max_workers=max_workers)
    for out, map_f in outputs:
        for idx, v in enumerate(out):
            out_f.write(
                f'{idx}\t{v}\t{map_f["algorithm"]}\t{map_f["randomSeed"]}\n')
    out_f.close()
    dfs = []
    for out, map_f in outputs2:
        out['algorithm'] = map_f["algorithm"]
        out['randomSeed'] = map_f["randomSeed"]
        dfs.append(out)
    dfs = pd.concat(dfs)
    dfs.to_csv(out_eci, sep='\t', header=None, index=None)
