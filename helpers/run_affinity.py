# Compute group affinity and homophily in directed hypergraphs.
import sys
sys.path.insert(1, '../')
from helpers.config import max_workers, max_size, min_size
from helpers.config import sample_path, met_path, folder
from helpers.utils import compute_v_hlist_maps
from helpers.io import read_as_directed_hyperedge_list, read_classes
from helpers.parser_for_output_file_names import parser_for_output_file_name as parse_file_name
import numpy as np
import zipfile
import io
from collections import defaultdict
from tqdm.contrib.concurrent import process_map
from tqdm import tqdm


def affinity_head_one(heads,
                      tails,
                      classes,
                      num_classes,
                      sampler=None,
                      ss=None):
    '''
    Method used to compute the group affinity
    when the head size of each hyperedge is 1.

    INPUT
    =====
    heads, list: list of hyperedge heads.
    tails, list: list of hyperedge tails.
    classes, dict: for each vertex, its label.
    num_classes, int: number of distinct classes.
    sampler, str: name of the sampler used to generate the input hypergraph.
                  If None, the input is the observed hypergraph.
    ss, int: seed used by the sampler to generate the input hypergraph.
             If None, the input is the observed hypergraph.
    '''

    _, v_to_tlist, hedge_size_map = compute_v_hlist_maps(heads, tails)
    max_s = max([1 + len(tails[i]) for i in range(len(tails))])
    max_s = min(max_s, max_size)

    output = []

    for s in tqdm(range(min_size, max_s)):

        hyperedge_to_class_vs = np.zeros((len(heads), num_classes, 2),
                                         dtype=np.int32)
        for idx, h in enumerate(heads):
            if len(h) + len(tails[idx]) == s:
                for v in h:
                    hyperedge_to_class_vs[idx][classes[v]][0] += 1
                for v in tails[idx]:
                    hyperedge_to_class_vs[idx][classes[v]][1] += 1

        # score 1: probability that v is X given that head is X
        n_classes = np.zeros(num_classes)
        d_classes = np.zeros(num_classes)

        for v in v_to_tlist:
            x = classes[v]
            for h in v_to_tlist[v]:
                if hedge_size_map[h] != s:
                    continue
                if hyperedge_to_class_vs[h][x][0] == 1:
                    n_classes[x] += 1
                d_classes[x] += 1

        head_affinity_scores = np.divide(n_classes, d_classes)
        for idx, aff in enumerate(head_affinity_scores):
            # sampler, seed, class, size, beta, head/tail, affinity score
            if sampler is not None:
                output.append([sampler, ss, idx, s, 1, 'head', aff])
            else:
                output.append([idx, s, 1, 'head', aff])

        for beta in range(1, s + 1):

            # score 2: prob that v is X given that beta nodes are X in tail
            n_classes = np.zeros(num_classes)
            d_classes = np.zeros(num_classes)

            for v in v_to_tlist:
                x = classes[v]
                for h in v_to_tlist[v]:
                    if hedge_size_map[h] != s:
                        continue
                    if hyperedge_to_class_vs[h][x][1] == beta:
                        n_classes[x] += 1
                    d_classes[x] += 1

            tail_affinity_scores = np.divide(n_classes, d_classes)
            for idx, aff in enumerate(tail_affinity_scores):
                if sampler is not None:
                    output.append([sampler, ss, idx, s, beta, 'tail', aff])
                else:
                    output.append([idx, s, beta, 'tail', aff])
    return output


def homophily_head_one_all_classes(heads,
                                   tails,
                                   classes,
                                   num_classes,
                                   sampler=None,
                                   ss=None):
    '''
    Multi-class method used to compute the homophily 
    when the head size of each hyperedge is 1.

    INPUT
    =====
    heads, list: list of hyperedge heads.
    tails, list: list of hyperedge tails.
    classes, dict: for each vertex, its label.
    num_classes, int: number of distinct classes.
    sampler, str: name of the sampler used to generate the input hypergraph.
                  If None, the input is the observed hypergraph.
    ss, int: seed used by the sampler to generate the input hypergraph.
             If None, the input is the observed hypergraph.
    '''

    output = []

    HM = np.zeros((num_classes, num_classes), dtype=np.float64)

    for idx, h in enumerate(heads):
        c = classes[h[0]]
        for v in tails[idx]:
            HM[c][classes[v]] += (1. / len(tails[idx]))
    for i in range(num_classes):
        for j in range(num_classes):
            if sampler is not None:
                output.append([sampler, ss, f'M{i}-{j}', HM[i][j]])
            else:
                output.append([f'M{i}-{j}', HM[i][j]])
    return output


def homophily_head_one(heads,
                       tails,
                       classes,
                       sampler=None,
                       ss=None):
    '''
    Binary-class method used to compute the homophily 
    when the head size of each hyperedge is 1.

    INPUT
    =====
    heads, list: list of hyperedge heads.
    tails, list: list of hyperedge tails.
    classes, dict: for each vertex, its label.
    num_classes, int: number of distinct classes.
    sampler, str: name of the sampler used to generate the input hypergraph.
                  If None, the input is the observed hypergraph.
    ss, int: seed used by the sampler to generate the input hypergraph.
             If None, the input is the observed hypergraph.
    '''

    output = []
    counts_list = []

    # for each hyperedge, count the number
    # of nodes for each class
    for idx, h in enumerate(heads):
        #  we consider only Rep and Dem
        counts = np.zeros(2)
        c = classes[h[0]]
        if c not in [0, 1]:
            continue
        for v in tails[idx]:
            c2 = classes[v]
            if c2 not in [0, 1]:
                continue
            counts[c2] += 1.
        # class of the head, array of counts
        # print(c, counts)
        counts_list.append((c, counts))
    # sum over classes
    HM = np.zeros((2, 2), dtype=np.float64)
    for c1, count in counts_list:
        for c2 in range(2):
            if (count[0] + count[1]) > 0:
                HM[c1][c2] += count[c2] / (count[0] + count[1])
    # output
    for i in range(2):
        for j in range(2):
            if sampler is not None:
                output.append([sampler, ss, f'M{i}-{j}', HM[i][j]])
            else:
                output.append([f'M{i}-{j}', HM[i][j]])
    return output


def affinity_head_beta(heads,
                       tails,
                       classes,
                       num_classes,
                       sampler=None,
                       ss=None):
    '''
    Method used when the heads have heterogeneous sizes >= 1.

    INPUT
    =====
    heads, list: list of hyperedge heads.
    tails, list: list of hyperedge tails.
    classes, dict: for each vertex, its label.
    num_classes, int: number of distinct classes.
    sampler, str: name of the sampler used to generate the input hypergraph.
                  If None, the input is the observed hypergraph.
    ss, int: seed used by the sampler to generate the input hypergraph.
             If None, the input is the observed hypergraph.
    '''
    v_to_hlist, v_to_tlist, _ = compute_v_hlist_maps(heads, tails)
    # max hyperedge size to consider
    max_s = max([len(heads[i]) + len(tails[i]) for i in range(len(heads))])
    max_s = min(max_s, max_size)
    # max head size to consider
    max_h = max([len(h) for h in heads])

    output = []

    for s in tqdm(range(min_size, max_s)):

        hyperedge_to_class_vs = np.zeros((len(heads), num_classes, 2),
                                         dtype=np.int32)
        for idx, h in enumerate(heads):
            if len(h) + len(tails[idx]) == s:
                for v in h:
                    hyperedge_to_class_vs[idx][classes[v]][0] += 1
                for v in tails[idx]:
                    hyperedge_to_class_vs[idx][classes[v]][1] += 1

        for h_size in range(1, min(max_h, s) + 1):

            # a elements of a certain class in the head
            for a in range(1, h_size + 1):

                # probability that v is part of a tail of a hyperedge
                # where the head has a nodes with the same class
                n_classes = np.zeros(num_classes)
                d_classes = np.zeros(num_classes)

                for v in v_to_tlist:
                    x = classes[v]
                    for h in v_to_tlist[v]:
                        if len(heads[h]) != h_size:
                            continue
                        if hyperedge_to_class_vs[h][x][0] == a:
                            n_classes[x] += 1
                        d_classes[x] += 1
                head_affinity_scores = np.divide(n_classes, d_classes)
                for idx, aff in enumerate(head_affinity_scores):
                    if sampler is not None:
                        output.append([sampler, ss, idx,
                                       s, a, h_size, 'head', aff])
                    else:
                        output.append([idx, s, a, h_size, 'head', aff])

            # a elements of a certain class in the tail
            for a in range(1, s - h_size + 1):

                # probability that v is part of a head of a hyperedge
                # where the tail has a nodes with the same class
                n_classes = np.zeros(num_classes)
                d_classes = np.zeros(num_classes)

                for v in v_to_hlist:
                    x = classes[v]
                    for h in v_to_hlist[v]:
                        if len(heads[h]) != h_size:
                            continue
                        if hyperedge_to_class_vs[h][x][1] == a:
                            n_classes[x] += 1
                        d_classes[x] += 1
                tail_affinity_scores = np.divide(n_classes, d_classes)
                for idx, aff in enumerate(tail_affinity_scores):
                    if sampler is not None:
                        output.append([sampler, ss, idx,
                                       s, a, h_size, 'tail', aff])
                    else:
                        output.append([idx, s, a, h_size, 'tail', aff])
    return output


def parallel_affinity(inp):
    '''
    Run affinity experiment for the sample read from file_name.
    '''
    file_path = inp[0]
    file_name = inp[1]
    vmap = inp[2]
    classes = inp[3]
    num_classes = inp[4]
    typ = inp[5]

    z = zipfile.ZipFile(file_path, "r")
    with z.open(file_name, 'r') as f1:
        map_fields = parse_file_name(file_name, typ)
        fobj = io.TextIOWrapper(f1, encoding='utf-8', newline='')
        heads, tails = read_as_directed_hyperedge_list(fobj, vmap)

    return affinity_sample(heads,
                           tails,
                           classes,
                           num_classes,
                           map_fields['algorithm'],
                           map_fields['randomSeed'])


def affinity_sample(heads,
                    tails,
                    classes,
                    num_classes,
                    sampler=None,
                    ss=None):
    '''
    If all the heads have size 1, the method calls affinity_head_one;
    Otherwise it computes the traditional affinity measure.

    INPUT
    ------
    heads, list: list of hyperedge heads.
    tails, list: list of hyperedge tails.
    classes, dict: for each vertex, its label.
    num_classes, int: number of distinct classes.
    sampler, str: name of the sampler used to generate the input hypergraph.
                  If None, the input is the observed hypergraph.
    ss, int: seed used by the sampler to generate the input hypergraph.
             If None, the input is the observed hypergraph.
    '''
    out2 = None
    if max([len(h) for h in heads]) == 1:
        out1 = affinity_head_one(heads, tails,
                                 classes, num_classes,
                                 sampler,
                                 ss)
        out2 = homophily_head_one(heads, tails,
                                  classes,
                                  sampler,
                                  ss)
    else:
        out1 = affinity_head_beta(heads, tails,
                                  classes, num_classes,
                                  sampler,
                                  ss)
    return (out1, out2)


if __name__ == '__main__':

    root = sys.argv[1]
    typ = sys.argv[2]
    data_path = folder + root
    print(root, typ)

    # run affinity observed hypergraph
    vmap = defaultdict(int)
    with open(f'{data_path}.tsv') as in_f:
        heads, tails = read_as_directed_hyperedge_list(in_f, vmap)

    # case: congress bills
    if "congress_bills" in root:
        class_path = folder + "congress_bills_classes.tsv"
    else:
        # case: dblp
        class_path = f'{data_path}_classes.tsv'
    # read node labels
    with open(class_path) as in_f:
        classes = read_classes(in_f, vmap)
    num_classes = max(classes.values()) + 1

    out1, out2 = affinity_sample(heads,
                                 tails,
                                 classes,
                                 num_classes)
    header3 = 'Type\tValue\n'
    if len(out1[0]) == 5:
        # case: head sizes = 1
        header = 'Class\tHyperedge Size\tAlpha\tType\tAffinity\n'
    else:
        # case: head sizes >= 1
        header = 'Class\tHyperedge Size\tAlpha\tBeta\tType\tAffinity\n'

    with open(met_path + f"affinity_{root}.tsv", 'w') as out_f:
        out_f.write(header)
        for o in out1:
            out_f.write('\t'.join([str(v) for v in o]) + '\n')
    if out2 is not None:
        with open(met_path + f"homophily_{root}.tsv", 'w') as out_f:
            out_f.write(header3)
            for o in out2:
                out_f.write('\t'.join([str(v) for v in o]) + '\n')

    file_path = f'{sample_path}/{typ}/{root}.zip'
    out_f = open(met_path + f"affinity_{root}_samples_{typ}.tsv", 'w')
    out_f2 = open(met_path + f"homophily_{root}_samples_{typ}.tsv", 'w')

    z = zipfile.ZipFile(file_path, "r")
    zinfo = z.namelist()
    inputs = []
    for file_name in zinfo:
        if file_name.startswith(".") or file_name.startswith("__MACOSX"):
            continue
        if not file_name.endswith('.tsv'):
            continue
        inputs.append([file_path,
                       file_name,
                       vmap,
                       classes,
                       num_classes,
                       typ])

    outputs = process_map(parallel_affinity,
                          inputs,
                          max_workers=max_workers)

    header2 = f'Sampler\tSample Id\t{header}'
    header4 = f'Sampler\tSample Id\t{header3}'
    out_f.write(header2)
    out_f2.write(header4)
    for out1, out2 in outputs:
        for o in out1:
            out_f.write('\t'.join([str(v) for v in o]) + '\n')
        if out2 is not None:
            for o in out2:
                out_f2.write('\t'.join([str(v) for v in o]) + '\n')
    out_f.close()
    out_f2.close()
