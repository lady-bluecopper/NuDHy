# Generate num_samples samples using the variant of ReDi
# given in input (values abdissible ReDi, ReDiD, Base, BaseD, Null)
import sys
sys.path.insert(1, '../')
from helpers.hyprec.HyperRec import measure_reciprocity
from helpers.hyprec.ReDi import parallel_redi, parallel_null
from helpers.io import read_as_directed_hyperedge_list
from helpers.config import sample_path, max_workers, alpha, num_samples, folder
import numpy as np
import time
from collections import defaultdict
from tqdm.contrib.concurrent import process_map


if __name__ == '__main__':
    root = sys.argv[1]
    typ = sys.argv[2]
    data_path = folder + root
    print(root, typ)

    # READ ORIGINAL
    vmap = defaultdict(int)
    with open(f'{data_path}.tsv') as in_f:
        heads, tails = read_as_directed_hyperedge_list(in_f, vmap)
    head_i = np.array([set(t) for t in heads])
    tail_i = np.array([set(t) for t in tails])
    # inv map
    inv_map = {v: k for k, v in vmap.items()}
    print(data_path)
    beta1 = 0
    beta2 = 0
    if (typ == 'ReDi') or (typ == 'ReDiD'):
        print('starting computation of reciprocity (ns).', time.time_ns())
        # reciprocity in observed dataset
        recip_vector = measure_reciprocity(head_index=head_i,
                                           tail_index=tail_i,
                                           alpha=alpha,
                                           special_case="normal")
        print('reciprocity observed dataset computed (ns).', time.time_ns())
        # compute beta1
        actually_recip = []
        for v in recip_vector:
            if v > 0:
                actually_recip.append(v)
        beta1 = float(len(actually_recip)) / len(head_i)
        # compute beta2
        beta2 = sum(actually_recip)
        if len(actually_recip) > 0:
            beta2 /= len(actually_recip)

    # GENERATE SAMPLES
    folder = f'{sample_path}/{typ}/{root}'
    # bashCommand = f'mkdir -p {folder}'
    # os.system(bashCommand)

    if typ == 'Null':
        inputs = [[head_i, tail_i,
                   i, folder, inv_map] for i in range(num_samples)]
        outputs = process_map(parallel_null, inputs, max_workers=max_workers)
    else:
        inputs = [[head_i, tail_i,
                   beta1, beta2, i,
                   folder, inv_map, typ] for i in range(num_samples)]
        outputs = process_map(parallel_redi, inputs, max_workers=max_workers)
