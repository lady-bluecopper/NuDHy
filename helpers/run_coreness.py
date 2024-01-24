import sys
sys.path.insert(1, '../')
from helpers.config import sample_path, max_workers, folder
from helpers.io import read_as_directed_hyperedge_list
from helpers.coreness import run_coreness
import numpy as np
import zipfile
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
from collections import defaultdict


if __name__ == '__main__':

    root = sys.argv[1]
    typ = sys.argv[2]
    data_path = folder + root
    print(root, typ)

    # READ ORIGINAL
    vmap = defaultdict(int)
    with open(f'{data_path}.tsv') as in_f:
        heads, tails = read_as_directed_hyperedge_list(in_f, vmap)
    head_index = np.array([set(t) for t in heads])
    tail_index = np.array([set(t) for t in tails])

    with open(f'{data_path}_vmap.tsv', 'w') as o_f:
        for k, v in vmap.items():
            o_f.write(f'{v}\t{k}\n')

    num_vertices = len(vmap)

    max_m = max([len(heads[i]) + len(tails[i]) for i in range(len(tails))]) + 1
    max_m = min(max_m, 10)

    # hyper-coreness in observed hypergraph
    run_coreness([data_path + '.tsv', root, vmap,
                  max_m, num_vertices, root, None])

    print('m-h-coreness and m-t-coreness observed hypergraph computed.')

    file_path = f'{sample_path}/{typ}/{root}.zip'

    z = zipfile.ZipFile(file_path, "r")
    zinfo = z.namelist()
    inputs = []
    for file_name in tqdm(zinfo):
        if file_name.startswith(".") or file_name.startswith("__MACOSX"):
            continue
        if not file_name.endswith('.tsv'):
            continue
        inputs.append([file_path, file_name, vmap,
                       max_m, num_vertices, root, typ])

    process_map(run_coreness, inputs, max_workers=max_workers)
    print('m-h-coreness and m-t-coreness samples computed.')
