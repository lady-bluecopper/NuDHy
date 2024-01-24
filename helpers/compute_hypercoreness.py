# Given the output files of run_coreness.py,
# it computes the hyper-coreness of each sample generated
# by the sampler given in input
import sys
sys.path.insert(1, '../')
from helpers.io import read_as_directed_hyperedge_list
from helpers.config import num_samples, met_path, folder
import numpy as np
from collections import defaultdict


def compute_hyper_coreness(shell_indices, num_v):
    '''
    INPUT
    ======
    shell_indices, dict: for each m, the array of m-shell values of each v
    num_v, int: number of vertices in the hypergraph

    OUTPUT
    ======
    hypercoreness, np.array: for each vertex, the hyper-coreness
    '''
    hypercoreness = np.zeros(num_v)
    max_coreness = dict()
    for m in shell_indices:
        max_coreness[m] = np.max(shell_indices[m])
    # compute hyper-coreness
    for m in shell_indices:
        if max_coreness[m] > 0:
            # normalize coreness values by max
            hypercoreness += np.array([1.0 * x / max_coreness[m]
                                      for x in shell_indices[m]])
    return hypercoreness


def get_shell_indices(file_name, num_v):
    with open(file_name) as in_f:
        # for each m, the array of m-shell values of each v
        shell_indices = dict()
        # for each m, for each k, the number of nodes
        # with m-shell index equal to k
        core_sizes = dict()
        sh_m = np.zeros(num_v)
        c_size = defaultdict(int)
        m = 2
        in_f.readline()
        for line in in_f.readlines():
            lst = line.split('\t')
            vidx = int(lst[0].strip())
            this_m = int(lst[1].strip())
            this_c = int(lst[2].strip())
            if this_m != m:
                # save shell-indices
                shell_indices[m] = sh_m
                sh_m = np.zeros(num_v)
                # save core sizes
                core_sizes[m] = c_size
                c_size = defaultdict(int)
                m = this_m
            sh_m[vidx] = this_c
            c_size[this_c] += 1
    return shell_indices, core_sizes


if __name__ == '__main__':

    root = sys.argv[1]
    typ = sys.argv[2]
    if typ.lower() == 'nudhy':
        samplers = ['NuDHy_A', 'NuDHy_C']
    else:
        samplers = [typ]

    data_path = folder + root
    print(root, typ)

    # READ ORIGINAL
    vmap = defaultdict(int)
    with open(f'{data_path}.tsv') as in_f:
        heads, tails = read_as_directed_hyperedge_list(in_f, vmap)
    num_vertices = len(vmap)

    for oo in ['o', 'i']:
        # OBSERVED HYPERGRAPH
        core_path = f'{met_path}{root}_coreness_{oo}.tsv'
        core_s_path = f'{met_path}{root}_core_size_{oo}.tsv'
        hcore_path = f'{met_path}{root}_hypercoreness_{oo}.tsv'
        shell_indices, core_sizes = get_shell_indices(core_path,
                                                      num_vertices)
        # SAVE CORE SIZES
        with open(core_s_path, 'w') as out_f:
            out_f.write('m\tm-shellIndex\tSize\n')
            for m in core_sizes:
                for k, size in core_sizes[m].items():
                    out_f.write(f'{m}\t{k}\t{size}\n')
        # COMPUTE HYPERCORENESS FOR OBSERVED HYPERGRAPH
        hypercoreness = compute_hyper_coreness(shell_indices,
                                               num_vertices)
        # SAVE HYPER-CORENESS
        with open(hcore_path, 'w') as out_f:
            out_f.write('VertexId\tHyper-coreness\n')
            for idx, hc in enumerate(hypercoreness):
                out_f.write(f'{idx}\t{hc}\n')
        # SAMPLES
        for sampler in samplers:
            for sample in np.arange(num_samples):

                core_path = f'{met_path}{root}__randomSeed_{sample}'
                core_path += f'__algorithm_{sampler}__coreness_{oo}.tsv'

                core_s_path = f'{met_path}{root}__randomSeed_{sample}'
                core_s_path += f'__algorithm_{sampler}__core_size_{oo}.tsv'

                hcore_path = f'{met_path}{root}__randomSeed_{sample}'
                hcore_path += f'__algorithm_{sampler}__hypercoreness_{oo}.tsv'

                shell_indices, core_sizes = get_shell_indices(
                    core_path, num_vertices)

                # SAVE CORE SIZES FOR THIS SAMPLE
                with open(core_s_path, 'w') as out_f:
                    out_f.write('m\tm-shellIndex\tSize\n')
                    for m in core_sizes:
                        for k, size in core_sizes[m].items():
                            out_f.write(f'{m}\t{k}\t{size}\n')
                # HYPERCORENESS FOR SAMPLES
                hypercoreness = compute_hyper_coreness(
                    shell_indices, num_vertices)
                with open(hcore_path, 'w') as out_f:
                    out_f.write('VertexId\tHyper-coreness\n')
                    for idx, hc in enumerate(hypercoreness):
                        out_f.write(f'{idx}\t{hc}\n')
