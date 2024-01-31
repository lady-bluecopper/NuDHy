# Compute GENEPY and ECI for directed bipartite graphs representing
# directed hypergraphs
import numpy as np
import econci
import pandas as pd
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import eigsh
from copy import deepcopy
from tqdm.contrib.concurrent import process_map
import zipfile
import io
from datetime import datetime
import sys
sys.path.insert(1, '../')
from helpers.config import folder, side, group, met_path, sample_path, max_workers
from helpers.io import read_as_bipartite
from helpers.parser_for_output_file_names import parser_for_output_file_name as parse_file_name


def fitness_complexity_calculation(m, iterations=1000, samples=0):
    """
    @Author
    Orazio Angelini
    @Paper
    Angelini, O., Cristelli, M., Zaccaria, A., & Pietronero, L. (2017).
    The complex dynamics of products and its asymptotic properties.
    PLoS ONE, 12(5).
    @GitHub
    https://github.com/ganileni/ectools

    Returns fitness and complexity at different iterations for a matrix m.
    m is intended to be of shape [n_countries, n_products], i.e.
    Fitness will be computed on the rows and Complexity on the columns of m.
    Note: this function handles all data cleaning procedures.
    The actual fixed point calculation is done by the function iterative_map().

    Arguments:
        m, np.ndarray: binary matrix of shape [n_countries, n_products]
        iterations, int: iterations to perform
        samples, int: samples of fitness and complexity to take while iterating

    Returns:
        fit_v, np.ndarray: matrix [samples, n_countries] with the
                values of Fitness for each country taken while iterating
        fit_norm_v, np.ndarray: array [samples] with the
                normalized Fitness measured while iterating
        compl_v, np.ndarray: matrix [samples, n_products] with the
                values of Complexity for each product taken while iterating
        compl_norm_v, np.ndarray: array [samples] with the
                normalized Complexity measured while iterating
        sample_times, np.ndarray: array [samples] containing the iteration at
                which the Fitness and Complexity values were measured

    """
    assert np.isnan(m).sum() == 0
    # if there are some products exported by nobody (!)
    noones_products = m.sum(axis=0) == 0
    noones_products_present = noones_products.sum()
    if noones_products_present:
        # remove columns corresponding to these products
        m = m[:, np.invert(noones_products)]
    # if there are some countries that export nothing (!)
    nontrade_countries = m.sum(axis=1) == 0
    nontrade_countries_present = nontrade_countries.sum()
    if nontrade_countries_present:
        # remove rows corresponding to such countries
        m = m[np.invert(nontrade_countries), :]
    # how big is the matrix
    num_countries, num_products = m.shape
    # where to take samples:
    if samples != 0:
        step = int(iterations / samples)
        sample_times = np.zeros(samples)
    else:
        step = iterations
        sample_times = np.array([0])
    # initialize fitness and complexity
    fitness = np.array(np.ones(num_countries), dtype=np.float64)
    complexity = np.array(np.ones(num_products), dtype=np.float64)
    # initialize sampling lists so that the algorithm
    # doesn't have to do dynamic resizing of arrays
    fit_v = np.zeros([samples + 1, num_countries])
    fit_norm_v = np.zeros(samples + 1)
    compl_v = np.zeros([samples + 1, num_products])
    compl_norm_v = np.zeros(samples + 1)
    fit_v, fit_norm_v, compl_v, compl_norm_v, sample_times = \
        iterative_map(m, iterations, step,
                      fitness, complexity,
                      fit_v, fit_norm_v,
                      compl_v, compl_norm_v,
                      sample_times)
    # values for countries that export nothing reinserted as nans:
    # (must come before reinsertion of problematic complexity values!)
    if nontrade_countries_present:
        # make empty arrays with adequate size
        fit_v_complete = np.zeros([
            fit_v.shape[0],
            len(nontrade_countries)])
        nan_vector = np.zeros(fit_v_complete.shape[0])
        nan_vector[:] = np.nan
        # iteratively reinsert nontrading countries
        # preserving original ordering
        count = 0
        for k, nontrading in enumerate(nontrade_countries):
            # if it's a trading country, put in its fitness
            if not nontrading:
                fit_v_complete[:, k] = fit_v[:, count]
                count += 1
            # otherwise put in a nan column
            else:
                fit_v_complete[:, k] = nan_vector
        # substitute complete values in vector to be returned
        fit_v = fit_v_complete
    # values for products exported by nobody reinserted as nans:
    if noones_products_present:
        # make empty arrays with adequate size
        compl_v_complete = np.zeros([
            compl_v.shape[0],
            len(noones_products)])
        nan_vector = np.zeros([compl_v_complete.shape[0]])
        nan_vector[:] = np.nan
        # iteratively reinsert noones products preserving original ordering
        count = 0
        for k, noones in enumerate(noones_products):
            # if it's not a problematic product
            if not noones:
                compl_v_complete[:, k] = compl_v[:, count]
                count += 1
            else:
                compl_v_complete[:, k] = nan_vector
        # substitute complete values in vector to be returned
        compl_v = compl_v_complete
    return fit_v, fit_norm_v, compl_v, compl_norm_v, sample_times


def iterative_map(m, iterations, step,
                  fitness, complexity,
                  fit_v,
                  fit_norm_v,
                  compl_v,
                  compl_norm_v, sample_times):
    '''
    m: Binary matrix on which to compute the Fitness-Complexity algorithm.
    '''
    """
    @Author
    Orazio Angelini
    @Paper
    Angelini, O., Cristelli, M., Zaccaria, A., & Pietronero, L. (2017).
    The complex dynamics of products and its asymptotic properties.
    PLoS ONE, 12(5).
    @GitHub
    https://github.com/ganileni/ectools

    Arguments:
        m, np.ndarray: binary matrix of shape [n_countries, n_products]
        iterations, int: iterations to perform
        step, int: step at which fitness and complexity must be output
        fitness, np.ndarray: array [n_countries] where the fitness values
                will be iteratively updated
        complexity, np.ndarray: array [n_products] where the complexity values
                will be iteratively updated
        samples, int: samples of fitness and complexity to take while iterating
        fit_v, np.ndarray: matrix [samples, n_countries] to store the
                values of Fitness for each country taken after step iterations
        fit_norm_v, np.ndarray: array [samples] to store the values of
                normalized Fitness measured after step iterations
        compl_v, np.ndarray: matrix [samples, n_products] to store the values
                of Complexity for each product taken after step iterations
        compl_norm_v, np.ndarray: array [samples] to store the values of
                normalized Complexity measured after step iterations
        sample_times, np.ndarray: array [samples] to store the iteration at
                which the Fitness and Complexity values were measured
    Returns:
        fit_v, np.ndarray: matrix [samples, n_countries] with the
                values of Fitness for each country taken while iterating
        fit_norm_v, np.ndarray: array [samples] with the
                normalized Fitness measured while iterating
        compl_v, np.ndarray: matrix [samples, n_products] with the
                values of Complexity for each product taken while iterating
        compl_norm_v, np.ndarray: array [samples] with the
                normalized Complexity measured while iterating
        sample_times, np.ndarray: array [samples] containing the iteration at
                which the Fitness and Complexity values were measured

    """
    samples_taken = 0
    for n in range(iterations):
        tmp_fitness = fitness
        tmp_complexity = complexity
        fitness = np.dot(m, tmp_complexity)
        complexity = 1 / (np.dot(1 / tmp_fitness, m))
        # normalize
        fitness = fitness / (fitness.sum() / len(fitness))
        complexity = complexity / (complexity.sum() / len(complexity))
        # check for nans in complexity
        # (nans appear there first, then in fitness)
        if np.isnan(complexity).sum():
            # remove all subsequent values from samples arrays
            fit_v = fit_v[:samples_taken + 1]
            compl_v = compl_v[:samples_taken + 1]
            fit_norm_v = fit_norm_v[:samples_taken + 1]
            compl_norm_v = compl_norm_v[:samples_taken + 1]
            sample_times = sample_times[:samples_taken + 1]
            # recover last non-nan fitness and complexity
            complexity = tmp_complexity
            fitness = tmp_fitness
            break
        # save iterations
        if n % step == 0:
            fit_v[samples_taken] = fitness
            compl_v[samples_taken] = complexity
            fit_norm_v[samples_taken] = fitness.sum()
            compl_norm_v[samples_taken] = complexity.sum()
            sample_times[samples_taken] = n
            samples_taken += 1
    fit_v[-1] = fitness
    compl_v[-1] = complexity
    fit_norm_v[-1] = fitness.sum()
    compl_norm_v[-1] = complexity.sum()
    sample_times[-1] = n
    return fit_v, fit_norm_v, compl_v, compl_norm_v, sample_times


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


def parallel_fitness(inp):
    '''
    Run Fitness for the sample read from file_name
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
        M = biadjacency_mat(fobj, vmap, side).toarray()
    fit_matrix, _, _, _, _ = fitness_complexity_calculation(M)
    return fit_matrix, map_fields


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
    print(datetime.now())
    vmap = dict()
    with open(f'{data_path}.tsv') as in_f:
        M = biadjacency_mat(in_f, vmap, side)
    _, _, genepy = genepy_index(M, group)
    with open(met_path + f"genepy_{root}_{side}.tsv", 'w') as out_f:
        for idx, v in enumerate(genepy):
            out_f.write(f'{idx}\t{v}\n')
    print('GENEPY in original hypergraph computed.')
    # run ECI observed hypergraph
    print(datetime.now())
    with open(f'{data_path}.tsv') as in_f:
        df = create_country_product_df(in_f, vmap, side)
    comp = econci.Complexity(df, c='country', p='product', values='export')
    comp.calculate_indexes()
    out_file = met_path + f"eci_{root}_{side}.tsv"
    comp.eci.to_csv(out_file, sep="\t", header=False, index=False)
    print('ECI in original hypergraph computed.')
    # run FITNESS observed hypergraph
    print(datetime.now())
    fit_map, _, _, _, _ = fitness_complexity_calculation(M.toarray())
    with open(met_path + f"fitness_{root}_{side}.tsv", 'w') as out_f:
        for idx, v in enumerate(fit_map[-1]):
            out_f.write(f'{idx}\t{v}\n')
    print('Fitness in original hypergraph computed.')
    # save index
    with open(met_path + f"vmap_{root}.tsv", 'w') as o_f:
        for key in vmap:
            o_f.write(f'{key}\t{vmap[key]}\n')
    # SAMPLES
    print(datetime.now())
    print('Starting analysis in samples...')
    file_path = f'{sample_path}/{typ}/{root}.zip'
    out_f = open(met_path + f"genepy_{root}_samples_{typ}_{side}.tsv", 'w')
    out_eci = met_path + f"eci_{root}_samples_{typ}_{side}.tsv"
    out_fit = open(met_path + f"fitness_{root}_samples_{typ}_{side}.tsv", 'w')

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
    # FITNESS
    outputs3 = process_map(parallel_fitness,
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
    for out, map_f in outputs3:
        for idx, v in enumerate(out[-1]):
            out_fit.write(
                f'{idx}\t{v}\t{map_f["algorithm"]}\t{map_f["randomSeed"]}\n')
    out_fit.close()
