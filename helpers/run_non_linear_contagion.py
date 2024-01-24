import sys
sys.path.insert(1, '../')
import time
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
from copy import deepcopy
import io
import zipfile
import numpy as np
from helpers.config import burnin_dec_dt, measure_dec_dt, nl_params
from helpers.config import lower_param, nb_points, recovery_rate, nb_history, init_I_list
from helpers.config import sample_path, met_path, folder, max_workers, seed
from helpers.io import process_data_nl_contagion
from helpers.parser_for_output_file_names import parser_for_output_file_name as parse_file_name
from _schon import PowerlawGroupSIS
from gcm import invasion_threshold_safe


# infection function
def beta(n, i, scale, shape):
    return scale * i**shape


def get_bounds(beta, nmax, scale, shape):
    max_inf = nmax * beta(nmax, nmax, scale, shape)  # upper bound
    min_inf = beta(1, 1, scale, shape)  # lower bound
    upper = max((max_inf, 1))
    lower = min((min_inf, 1))
    return (lower, upper)


def experience(inp):
    '''
    Simulate contagion process for the dataset in edge_list
    '''
    gm = inp[0]
    pn = inp[1]
    nmax = inp[2]
    edge_list = inp[3]
    shapes = inp[4]
    higher_param = inp[5]
    burnin_dt = inp[6]
    measure_dt = inp[7]
    quasi = inp[8]
    ss = inp[9]

    output = dict()

    for shape_infection in tqdm(shapes):
        print('shape infection', shape_infection)
        # get invasion threshold and bistability threshold
        scale_c = invasion_threshold_safe(beta, gm, pn,
                                          fixed_args=(shape_infection,),
                                          min_param=10**(-14), max_param=1)
        # define simulation parameters and output
        low_param = lower_param * scale_c
        high_param = higher_param * scale_c
        param_list = np.linspace(low_param, high_param, nb_points)
        results = dict()

        for init_I_upper in tqdm(init_I_list):
            print('init_I_upper', init_I_upper)
            results[init_I_upper] = dict()
            for scale in param_list:
                now = time.time()
                # varying parameters
                rate_bounds = get_bounds(beta, nmax, scale, shape_infection)
                I_list = []
                MI_list = []
                IS_list = []
                results[init_I_upper][scale] = dict()
                # define process
                cont = PowerlawGroupSIS(edge_list, recovery_rate, scale,
                                        shape_infection, rate_bounds)
                cont.infect_fraction(init_I_upper)
                cont.seed(ss)
                cont.initialize_history(nb_history)

                # define some measures
                cont.measure_prevalence()
                cont.measure_marginal_infection_probability()
                # cont.measure_infectious_set()

                # evolve in the quasistationary state
                # without measuring (burn-in)
                cont.evolve(burnin_dt, burnin_dec_dt, measure=False,
                            quasistationary=quasi)
                # evolve and measure
                cont.evolve(measure_dt, measure_dec_dt, measure=True,
                            quasistationary=quasi)

                # get the measure
                for measure in cont.get_measure_vector():
                    name = measure.get_name()
                    if name == "prevalence":
                        I_list = measure.get_result()
                    elif name == "marginal_infection_probability":
                        MI_list = measure.get_result()
                    # elif name == "infectious_set":
                    #     IS_list = measure.get_result()
                results[init_I_upper][scale]['prevalence'] = I_list
                results[init_I_upper][scale]['marginal_infection_probability'] = MI_list
                results[init_I_upper][scale]['infectious_set'] = IS_list
                print('examined', scale, 'in', (time.time() - now))
        output[shape_infection] = results
    return output


def parallel_simulations(inp):
    '''
    Pre-processes the sample in file_name and runs
    the simulation of the contagion process
    '''
    file_path = inp[0]
    file_name = inp[1]
    shapes = inp[2]
    higher_param = inp[3]
    burnin_dt = inp[4]
    measure_dt = inp[5]
    quasi = inp[6]
    ss = inp[7]
    typ = inp[8]

    z = zipfile.ZipFile(file_path, "r")
    with z.open(file_name, 'r') as f1:
        map_f = parse_file_name(file_name, typ)
        fobj = io.TextIOWrapper(f1, encoding='utf-8', newline='')
        output = process_data_nl_contagion(fobj)

    args = [output['gm'],
            output['pn'],
            output['nmax'],
            output['edge_list'],
            shapes,
            higher_param,
            burnin_dt,
            measure_dt,
            quasi,
            ss]
    return (experience(args), map_f['algorithm'], map_f['randomSeed'])


def write_results(res, seed, quasi, f1, f2, f3, ss=None, algo=None):
    # Write the results of the simulation process on disk
    if ss is not None:
        row_temp = [ss, algo]
    else:
        row_temp = []

    for sh_inf, out_dic in res.items():
        # initial density -> dict of results
        for I0_perc, out_dic_I0 in out_dic.items():
            # out_dic_I0 contains a key for each scale in param_list
            for sc in out_dic_I0:
                # 'prevalence' -> list
                # 'marginal_infection_probability' -> list
                # 'infectious_set' -> list
                I_list = out_dic_I0[sc]['prevalence']
                for v in I_list:
                    row = deepcopy(row_temp)
                    row.extend([sh_inf, I0_perc, sc, v, seed, quasi])
                    f1.write('\t'.join([str(e) for e in row]) + '\n')
                MI_list = out_dic_I0[sc]['marginal_infection_probability']
                for v in MI_list:
                    row = deepcopy(row_temp)
                    row.extend([sh_inf, I0_perc, sc, v, seed, quasi])
                    f2.write('\t'.join([str(e) for e in row]) + '\n')
                IS_list = out_dic_I0[sc]['infectious_set']
                for v in IS_list:
                    row = deepcopy(row_temp)
                    row.extend([sh_inf, I0_perc, sc, v, seed, quasi])
                    f3.write('\t'.join([str(e) for e in row]) + '\n')


if __name__ == '__main__':

    root = sys.argv[1]
    typ = sys.argv[2]
    data_path = folder + root
    print(root, typ)

    root_params = nl_params[root]
    shapes = root_params['shape_infection']
    higher_param = root_params['higher_param']
    burnin_dt = root_params['burnin_dt']
    measure_dt = root_params['measure_dt']
    quasi = root_params['quasi']

    header = 'Shape Infection\tInitial Density Infected\tLambda'
    header1 = header + '\tNum Infected\tSeed\tQuasi\n'
    header2 = header + '\tMarginal Infection Prob\tSeed\tQuasi\n'
    header3 = header + '\tInfectious Set\tSeed\tQuasi\n'
    header1s = 'Sample Id\tSampler\t' + header1
    header2s = 'Sample Id\tSampler\t' + header2
    header3s = 'Sample Id\tSampler\t' + header3

    # run contagion on observed hypergraph
    with open(f'{data_path}.tsv') as in_f:
        output = process_data_nl_contagion(in_f)

    args = [output['gm'],
            output['pn'],
            output['nmax'],
            output['edge_list'],
            shapes,
            higher_param,
            burnin_dt,
            measure_dt,
            quasi,
            seed]
    print('Starting simulations...')
    results = experience(args)
    out_f1 = open(met_path + f"rhos_{root}_dt{burnin_dt}_q{quasi}.tsv", 'w')
    out_f1.write(header1)
    out_f2 = open(met_path + f"I_prob_{root}_dt{burnin_dt}_q{quasi}.tsv", 'w')
    out_f2.write(header2)
    out_f3 = open(met_path + f"I_sets_{root}_dt{burnin_dt}_q{quasi}.tsv", 'w')
    out_f3.write(header3)

    write_results(results, seed, quasi, out_f1, out_f2, out_f3)

    out_f1.close()
    out_f2.close()
    out_f2.close()

    out_f1 = open(
        met_path + f"rhos_{root}_samples_{typ}_dt{burnin_dt}_q{quasi}.tsv", 'a')
    out_f1.write(header1s)
    out_f2 = open(
        met_path + f"I_prob_{root}_samples_{typ}_dt{burnin_dt}_q{quasi}.tsv", 'a')
    out_f2.write(header2s)
    out_f3 = open(
        met_path + f"I_sets_{root}_samples_{typ}_dt{burnin_dt}_q{quasi}.tsv", 'a')
    out_f3.write(header3s)

    file_path = f'{sample_path}/{typ}/{root}.zip'

    z = zipfile.ZipFile(file_path, "r")
    zinfo = z.namelist()
    inputs = []
    for file_name in zinfo:
        if file_name.startswith(".") or file_name.startswith("__MACOSX"):
            continue
        if not file_name.endswith('.tsv'):
            continue
        args = [file_path, file_name, shapes, higher_param,
                burnin_dt, measure_dt, quasi, seed, typ]
        inputs.append(args)

    all_results = process_map(parallel_simulations,
                              inputs,
                              max_workers=max_workers)
    for res, algo, ss in all_results:
        write_results(res, seed, quasi, out_f1, out_f2, out_f3, ss, algo)
    out_f1.close()
    out_f2.close()
    out_f2.close()
