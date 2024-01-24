# Methods to generate random directed hypergraphs using ReDi, ReDi-degreewise, and Null
from helpers.config import allow_inters
import numpy as np
from tqdm import tqdm
from itertools import combinations
import time
import sys
sys.path.insert(1, '../')


def naive_sampling(head_index, tail_index, seed):  # Naive Baseline Method!
    np.random.seed(seed)
    new_head = []
    new_tail = []

    head_size, tail_size = [], []

    v_max = 0

    total_head_nodes = 0
    total_tail_nodes = 0

    for i in range(head_index.shape[0]):

        cur_head_size = len(head_index[i])
        cur_tail_size = len(tail_index[i])

        total_head_nodes += cur_head_size
        total_tail_nodes += cur_tail_size

        head_size.append(cur_head_size)
        tail_size.append(cur_tail_size)
        union_max = max(head_index[i].union(tail_index[i]))
        if v_max < union_max:
            v_max = union_max

    in_degree = np.zeros(int(v_max + 1))
    out_degree = np.zeros(int(v_max + 1))

    for i in range(head_index.shape[0]):
        cur_head = list(head_index[i])
        cur_tail = list(tail_index[i])

        in_degree[cur_head] += 1
        out_degree[cur_tail] += 1

    total_degree = in_degree + out_degree
    required_index = np.where(total_degree != 0)[0]

    entire_nodes = np.arange(required_index.shape[0])

    entire_head_random = np.random.choice(entire_nodes, total_head_nodes)
    entire_tail_random = np.random.choice(entire_nodes, total_tail_nodes)

    head_indpt = 0
    tail_indpt = 0

    for i in (range(head_index.shape[0])):

        cond = True

        h_size = head_size[i]
        t_size = tail_size[i]

        next_head_indpt = head_indpt + h_size
        next_tail_indpt = tail_indpt + t_size

        cur_H = set(entire_head_random[head_indpt:next_head_indpt])
        cur_T = set(entire_tail_random[tail_indpt:next_tail_indpt])

        cond1 = len(cur_H) == h_size
        cond2 = len(cur_T) == t_size
        cond3 = len(cur_H.intersection(cur_T)) == 0

        if cond1 * cond2 * cond3 == 1:

            None

        else:
            while cond:

                cur_H = set(np.random.choice(
                    entire_nodes, h_size, replace=False))
                cur_T = set(np.random.choice(
                    entire_nodes, t_size, replace=False))

                if len(cur_H.intersection(cur_T)) == 0:
                    cond = False

        new_head.append((cur_H))
        new_tail.append((cur_T))
        head_indpt = next_head_indpt
        tail_indpt = next_tail_indpt

    return np.array(new_head), np.array(new_tail)


def find_interaction(total_tuple, valid_nodes, required_size, degree_list, delta):
    try:
        deg_lists = np.array(list(total_tuple[required_size].values()))
        probs = (deg_lists) / np.sum(deg_lists)
        c_inter = int(np.random.choice(a=np.arange(probs.shape[0]),
                                       size=1,
                                       p=probs,
                                       replace=False))
        add_list = np.array(list(total_tuple[required_size].keys()))[c_inter]
    except:
        valid_nodes = list(set(valid_nodes))
        deg_list = degree_list[:len(valid_nodes)].copy() + delta
        deg_p = deg_list / np.sum(deg_list)
        add_list = np.random.choice(a=valid_nodes,
                                    p=deg_p,
                                    size=required_size, replace=False)

    return set(add_list)


def update_interaction(total_tuple, new_edge):
    new_E = list(new_edge)

    for i in range(1, int(len(new_E) + 1)):
        # for each possible subset of i nodes
        # from new_E, increase its occurrence
        # count
        possible_comb = combinations(new_E, i)
        for t1 in possible_comb:
            try:
                total_tuple[i][t1] += 1
            except:
                total_tuple[i][t1] = 1
    return total_tuple


def get_recip_set(degree, candid_list, NPick, delta):
    deg_ = degree[list(candid_list)].copy()
    deg_ += delta
    prob = deg_ / np.sum(deg_)
    chosen_list = set(np.random.choice(a=list(candid_list),
                                       size=NPick,
                                       p=prob,
                                       replace=False))

    return chosen_list


def ReDi(head_index, tail_index,
         beta1, beta2,
         mode, seed=0,
         allow_inters=False):
    v_max = 0
    head_size = []
    tail_size = []
    np.random.seed(seed)
    delta = 1

    for i in range(head_index.shape[0]):
        # keep track of max node id
        union_max = max(head_index[i].union(tail_index[i]))
        if v_max < union_max:
            v_max = union_max
        if ((len(head_index[i])) <= 10) & ((len(tail_index[i])) <= 10):
            # store head and tail size of hyperedges with size <= 10
            head_size.extend([len(head_index[i])])
            tail_size.extend([len(tail_index[i])])
    # find node ids of nodes appearing in at least one hyperedge
    node_filter1 = np.zeros(int(v_max + 1))
    node_filter2 = np.zeros(int(v_max + 1))
    for heads, tails in zip(head_index, tail_index):
        node_filter1[list(heads)] += 1
        node_filter2[list(tails)] += 1
    condition_indicator = node_filter1 + node_filter2
    # compute number of such nodes
    total_N = np.where(condition_indicator != 0)[0].shape[0]
    # INITIALIZATION
    # these lists will store the heads and the tails
    # of the hyperedges of the sample
    total_heads = []
    total_tails = []
    # get max head and tail size
    max_S = max(np.max(head_size), np.max(tail_size))
    # these dictionaries will store, for each head and tail size i,
    # the number of occurrences of each subset of i nodes
    # encountered so far in total_heads and total_tails
    head_bucket = {i: {} for i in range(int(max_S + 1))}
    tail_bucket = {i: {} for i in range(int(max_S + 1))}
    # nodes added so far to the hypergraph
    current_nodes = []
    # stores the current in- and out-degrees of the nodes
    total_indegree = np.zeros(total_N)
    total_outdegree = np.zeros(total_N)
    # initialize graph with 2 * mode * max_S nodes
    init_n = max_S
    if mode > 1:
        init_n *= mode
    init_nodes = np.arange(int(init_n * 2))
    for i in range(int(init_n)):
        # initialize a new hyperedge with one node
        # in the head and a different node in the tail
        tmp_head_list = {init_nodes[int(i)]}
        tmp_tail_list = {init_nodes[int(init_n + i)]}
        # add this hyperedge to the hypergraph
        total_heads.append(tmp_head_list)
        total_tails.append(tmp_tail_list)
        # update head and tail subgroup counts
        head_bucket = update_interaction(total_tuple=head_bucket,
                                         new_edge=tmp_head_list)
        tail_bucket = update_interaction(total_tuple=tail_bucket,
                                         new_edge=tmp_tail_list)
        current_nodes.append(int(i))
        current_nodes.append(int(init_n + i))
        # update in- and out-degree of the nodes
        total_indegree[list(tmp_head_list)] += 1
        total_outdegree[list(tmp_tail_list)] += 1
    # INITIALIZATION BEFORE MAIN LOOP
    add_N = max(int(len(head_index) - total_N), 0)
    NP = np.ones(total_N)
    tmp_NP = np.ones(v_max + 1)

    if mode == 0:
        # select add_N random vertices uniformly at random
        add_pos = np.random.choice(a=np.arange(int(tmp_NP.shape[0])),
                                   size=add_N)
    else:
        # tmp_NP stores the number of hyperedges including nodes
        # with id lower than the id of this_node (used to find NP)
        for h_, t_ in zip(head_index, tail_index):
            this_node = max(max(h_), max(t_))
            tmp_NP[this_node] += 1
        # convert to probability
        tmpP = tmp_NP / np.sum(tmp_NP)
        np.random.shuffle(tmpP)
        # select add_N random vertices according to tmpP
        add_pos = np.random.choice(a=np.arange(int(tmpP.shape[0])),
                                   p=tmpP,
                                   size=add_N)
    # initialize vector that will tell us the number
    # of hyperedges that will contain each node
    for i in add_pos:
        NP[i] += 1
    each_add_N = NP
    # for each hyperedge, decide if it is going to be reciprocated
    # with probability beta1
    is_recip_vector = np.random.choice(a=[True, False],
                                       size=int(np.sum(each_add_N)),
                                       p=[beta1, 1 - beta1])
    # for each head, decide its size uniformly at random
    # from the head sizes of the hyperedges in the observed hypergraph
    size_vectors = np.random.choice(
        np.arange(len(head_size)), size=int(np.sum(each_add_N) * 3))
    # hyperedge id
    idx = -1
    # iterate over the nodes
    for i in tqdm(range(total_N)):
        # decide number of hyperedges that will include node i
        k = each_add_N[i]
        indptr = len(total_heads)
        # for each hyperedge to create
        for j in range(int(k)):

            cond = True
            idx += 1
            idx2 = idx - 1
            recip = is_recip_vector[idx]

            while cond:
                idx2 += 1
                cur_idx = size_vectors[idx2]

                N_head = int(head_size[cur_idx])
                N_tail = int(tail_size[cur_idx])
                # decide whether to insert the node in the
                # head or in the tail
                is_head = np.random.choice([True, False], 1)
                # if the hyperedge should be reciprocated
                if recip:
                    # If first iteration choose one edge and
                    # make it reciprocal to hyperedge j
                    if j == 0:
                        to_be_reciprocal = np.arange(len(total_heads))
                        recip_idx = int(np.random.choice(to_be_reciprocal, 1))
                        reciprocal_head = total_heads[recip_idx]
                        reciprocal_tail = total_tails[recip_idx]
                    else:
                        # candidate hyperedges among which we must choose
                        # the hyperedge recip_idx to make it reciprocal to hyperedge j
                        to_be_reciprocal = np.arange(indptr, len(total_heads))
                        recip_idx = int(np.random.choice(to_be_reciprocal, 1))
                        reciprocal_head = total_heads[recip_idx]
                        reciprocal_tail = total_tails[recip_idx]

                    max_tail_recips = min(N_tail, len(reciprocal_head)) - 1
                    max_head_recips = min(N_head, len(reciprocal_tail)) - 1

                    Nhead_reciprocal = 1
                    if max_head_recips > 1:
                        Nhead_reciprocal = np.random.binomial(n=max_head_recips,
                                                              p=beta2,
                                                              size=1) + 1
                    Ntail_reciprocal = 1
                    if max_tail_recips > 1:
                        Ntail_reciprocal = np.random.binomial(n=max_tail_recips,
                                                              p=beta2,
                                                              size=1) + 1

                    head_idx = set(get_recip_set(total_indegree,
                                                 reciprocal_tail,
                                                 Nhead_reciprocal,
                                                 delta))

                    tail_idx = set(get_recip_set(total_outdegree,
                                                 reciprocal_head,
                                                 Ntail_reciprocal,
                                                 delta))
                    if (len(head_idx) == N_head) and (len(tail_idx) != N_tail):
                        tail_idx = tail_idx.union({i})
                    elif len(head_idx) != N_head:
                        head_idx = head_idx.union({i})
                    # number of remaining nodes to add
                    # to the head and the tail of the hyperedges
                    N_add_head = N_head - len(head_idx)
                    N_add_tail = N_tail - len(tail_idx)
                    # fill the head if needed
                    if N_add_head > 0:
                        part_head = find_interaction(head_bucket,
                                                     current_nodes, N_add_head,
                                                     total_indegree, delta)
                        head_idx = head_idx.union(part_head)
                    # fill the tail if needed
                    if N_add_tail > 0:
                        part_tail = find_interaction(tail_bucket,
                                                     current_nodes, N_add_tail,
                                                     total_outdegree, delta)
                        tail_idx = tail_idx.union(part_tail)
                # this hyperedge should not have a reciprocal
                else:
                    # node i should be added to the head
                    if is_head:
                        part_head = find_interaction(head_bucket,
                                                     current_nodes, int(
                                                         N_head - 1),
                                                     total_indegree, delta)
                        head_idx = {i}.union(part_head)
                        tail_idx = find_interaction(tail_bucket, current_nodes, int(N_tail),
                                                    total_outdegree, delta)
                    # node i should be added to the tail
                    else:
                        head_idx = find_interaction(head_bucket, current_nodes, int(N_head),
                                                    total_indegree, delta)
                        part_tail = find_interaction(tail_bucket, current_nodes, int(N_tail - 1),
                                                     total_outdegree, delta)
                        tail_idx = {i}.union(part_tail)
                # the hyperedge is valid only if head and tail do not overlap
                if allow_inters or len(head_idx.intersection(tail_idx)) < 1:
                    cond = False
                    head_bucket = update_interaction(
                        head_bucket, head_idx)  # In-interaction
                    tail_bucket = update_interaction(
                        tail_bucket, tail_idx)  # Out-Interaction
                    total_indegree[list(head_idx)] += 1  # In-degree
                    total_outdegree[list(tail_idx)] += 1  # Out-degree
                    total_heads.append(head_idx)
                    total_tails.append(tail_idx)
        current_nodes.append(i)
    return np.array(total_heads), np.array(total_tails)


def ReDi_degreewise(head_index, tail_index,
                    beta1, beta2, mode,
                    seed=0, allow_inters=False):

    v_max = 0
    head_size = np.zeros(head_index.shape[0])
    tail_size = np.zeros(tail_index.shape[0])
    delta = 1

    np.random.seed(seed)

    for i in range(head_index.shape[0]):

        union_max = max(head_index[i].union(tail_index[i]))
        if v_max < union_max:
            v_max = union_max
        if len(head_index[i]) <= 10:
            head_size[i] = len(head_index[i])
        else:
            head_size[i] = 10
        if len(tail_index[i]) <= 10:
            tail_size[i] = len(tail_index[i])
        else:
            tail_size[i] = 10

    node_filter1 = np.zeros(int(v_max + 1))
    node_filter2 = np.zeros(int(v_max + 1))

    for heads, tails in zip(head_index, tail_index):
        node_filter1[list(heads)] += 1
        node_filter2[list(tails)] += 1

    condition_indicator = node_filter1 + node_filter2
    total_N = np.where(condition_indicator != 0)[0].shape[0]  # Number of nodes

    add_N = max(int(len(head_index) - total_N), 0)
    NP = np.ones(int(v_max) + 1)

    if mode == 0:  # mode 0 indicates naive addition of edges
        add_pos = np.random.choice(a=np.arange(int(total_N - 1)),
                                   size=add_N,
                                   replace=True)

        for i in add_pos:
            NP[i] += 1  # Add degree one by one

        each_add_N = np.random.choice(NP, total_N, replace=True)

    else:
        tmp_NP = np.ones(NP.shape[0])
        for h_, t_ in zip(head_index, tail_index):
            this_node = max(max(h_), max(t_))
            tmp_NP[this_node] += 1

        tmpP = tmp_NP / np.sum(tmp_NP)
        np.random.shuffle(tmpP)

        add_pos = np.random.choice(a=np.arange(int(tmpP.shape[0])),
                                   p=tmpP,
                                   size=add_N)

        for i in add_pos:
            NP[i] += 1  # Add degree one by one

        each_add_N = np.random.choice(NP, total_N, replace=True)

    total_heads = []
    total_tails = []

    in_norm = 0
    out_norm = 0

    # Initialization
    max_S = max(np.max(head_size), np.max(tail_size))
    if mode > 1:
        init_n = int(max_S) * mode
    else:
        init_n = int(max_S)
    init_nodes = np.arange(int(init_n * 2))

    current_nodes = []

    total_indegree = np.zeros(total_N)
    total_outdegree = np.zeros(total_N)
    IND = int(init_n * 2)

    for i in range(int(init_n)):
        tmp_head_list = {init_nodes[int(i)]}
        tmp_tail_list = {init_nodes[int(init_n + i)]}

        total_heads.append(tmp_head_list)
        total_tails.append(tmp_tail_list)

        in_norm += len(tmp_head_list)
        out_norm += len(tmp_head_list)
        total_indegree[list(tmp_head_list)] += 1
        total_outdegree[list(tmp_tail_list)] += 1

    is_recip_vector = np.random.choice(a=[True, False],
                                       size=int(np.sum(each_add_N)),
                                       p=[beta1, 1 - beta1])

    size_vectors = np.random.choice(np.arange(head_index.shape[0]),
                                    size=int(np.sum(each_add_N) * 3))
    idx = -1
    idx2 = -1

    for i in tqdm(range(each_add_N.shape[0])):

        if i > IND:
            IND = i

        k = each_add_N[i]  # Decide number of arcs to be created
        indptr = len(total_heads)
        for j in range(int(k)):
            idx += 1
            cond = True
            recip = is_recip_vector[idx]  # Check whether an arc is reciprocal
            idx2 = idx - 1
            while cond:
                idx2 += 1
                cur_idx = size_vectors[idx2]
                N_head = int(head_size[cur_idx])
                N_tail = int(tail_size[cur_idx])

                is_head = np.random.choice([True, False], 1)

                if recip:  # Let this arc reciprocal
                    # If first iteration ; Nothing but choice one edge and make as reciprocal

                    if j == 0:
                        to_be_reciprocal = np.arange(len(total_heads))
                        recip_idx = int(np.random.choice(to_be_reciprocal, 1))
                        reciprocal_head = total_heads[recip_idx]
                        reciprocal_tail = total_tails[recip_idx]

                    else:
                        to_be_reciprocal = np.arange(indptr, len(
                            total_heads))  # Reciprocal Candiate
                        # Choose one among them
                        recip_idx = int(np.random.choice(to_be_reciprocal, 1))
                        reciprocal_head = total_heads[recip_idx]
                        reciprocal_tail = total_tails[recip_idx]

                    max_tail_recips = min(N_tail, len(reciprocal_head)) - 1
                    max_head_recips = min(N_head, len(reciprocal_tail)) - 1

                    if max_head_recips > 1:
                        Nhead_reciprocal = np.random.binomial(n=max_head_recips,
                                                              p=beta2,
                                                              size=1) + 1
                    else:
                        Nhead_reciprocal = 1

                    if max_tail_recips > 1:
                        Ntail_reciprocal = np.random.binomial(n=max_tail_recips,
                                                              p=beta2,
                                                              size=1) + 1
                    else:
                        Ntail_reciprocal = 1

                    head_idx = set(get_recip_set(total_indegree,
                                                 reciprocal_tail,
                                                 Nhead_reciprocal,
                                                 delta))

                    tail_idx = set(get_recip_set(total_outdegree,
                                                 reciprocal_head,
                                                 Ntail_reciprocal,
                                                 delta))

                    if (len(head_idx) == N_head) & (len(tail_idx) == N_tail):

                        None  # Nothing to do

                    elif len(head_idx) == N_head:
                        tail_idx = tail_idx.union({i})

                    else:
                        head_idx = head_idx.union({i})

                    N_add_head = N_head - len(head_idx)
                    N_add_tail = N_tail - len(tail_idx)

                    if N_add_head > 0:
                        deg = total_indegree[:IND] + delta
                        prop = deg / np.sum(deg)
                        part_head = set(list(np.random.choice(a=np.arange(IND),
                                                              p=prop,
                                                              size=N_add_head,
                                                              replace=False)))
                        head_idx = head_idx.union(part_head)

                    if N_add_tail > 0:
                        deg = total_outdegree[:IND] + delta
                        prop = deg / np.sum(deg)
                        part_tail = set(list(np.random.choice(a=np.arange(IND),
                                                              p=prop,
                                                              size=N_add_tail,
                                                              replace=False)))
                        tail_idx = tail_idx.union(part_tail)

                else:  # Let this arc be naively created
                    if is_head:
                        deg = total_indegree[:IND] + delta
                        prop = deg / np.sum(deg)
                        part_head = set(list(np.random.choice(a=np.arange(IND),
                                                              p=prop,
                                                              size=int(
                                                                  N_head - 1),
                                                              replace=False)))
                        head_idx = {i}.union(part_head)

                        deg = total_outdegree[:IND] + delta
                        prop = deg / np.sum(deg)
                        tail_idx = set(list(np.random.choice(a=np.arange(IND),
                                                             p=prop,
                                                             size=int(N_tail),
                                                             replace=False)))
                    else:
                        deg = total_indegree[:IND] + delta
                        prop = deg / np.sum(deg)
                        head_idx = set(list(np.random.choice(a=np.arange(IND),
                                                             p=prop,
                                                             size=int(N_head),
                                                             replace=False)))

                        deg = total_outdegree[:IND] + delta
                        prop = deg / np.sum(deg)
                        part_tail = set(list(np.random.choice(a=np.arange(IND),
                                                              p=prop,
                                                              size=int(
                                                                  N_tail - 1),
                                                              replace=False)))
                        tail_idx = {i}.union(part_tail)

                if allow_inters or len(head_idx.intersection(tail_idx)) < 1:
                    cond = False
                    total_indegree[list(head_idx)] += 1  # In-degree
                    total_outdegree[list(tail_idx)] += 1  # Out-degree
                    total_heads.extend([head_idx])
                    total_tails.extend([tail_idx])

        current_nodes.append(i)

    return np.array(total_heads), np.array(total_tails)


def parallel_redi(inp):
    # head list of observed hypergraph
    headss = inp[0]
    # tail list of observed hypergraph
    tailss = inp[1]
    # proportion of reciprocal hyperedges
    beta1 = inp[2]
    # extent of reciprocity
    beta2 = inp[3]
    # seed for reproducibility
    seed = inp[4]
    # where to store the sample generated
    folder = inp[5]
    # inner-outer id map
    inv_map = inp[6]
    # type of sampler
    typ = inp[7]

    print('start generation (ns).', time.time_ns(), 'seed', seed)
    if (typ == 'ReDiD') or (typ == 'BaseD'):
        print('Running ReDi degree-wise')
        gen_head, gen_tail = ReDi_degreewise(head_index=headss,
                                             tail_index=tailss,
                                             beta1=beta1,
                                             beta2=beta2,
                                             mode=1,
                                             seed=seed,
                                             allow_inters=allow_inters)
    else:
        print('Running ReDi')
        gen_head, gen_tail = ReDi(head_index=headss,
                                  tail_index=tailss,
                                  beta1=beta1,
                                  beta2=beta2,
                                  mode=1,
                                  seed=seed,
                                  allow_inters=allow_inters)
    print('sample generated (ns).', time.time_ns(), 'seed', seed)

    with open(f"{folder}_{typ}{seed}_beta1{beta1}_beta2{beta2}.tsv", 'w') as out_f:
        if len(inv_map) > 0:
            for idx, h in enumerate(gen_head):
                hs = ','.join([str(inv_map[v]) for v in h])
                ts = ','.join([str(inv_map[v]) for v in gen_tail[idx]])
                out_f.write(f'{hs}\t{ts}\n')
        else:
            for idx, h in enumerate(gen_head):
                hs = ','.join([str(v) for v in h])
                ts = ','.join([str(v) for v in gen_tail[idx]])
                out_f.write(f'{hs}\t{ts}\n')


def parallel_null(inp):
    # head list of observed hypergraph
    headss = inp[0]
    # tail list of observed hypergraph
    tailss = inp[1]
    # seed for reproducibility
    seed = inp[2]
    # where to store the sample generated
    folder = inp[3]
    # inner-outer id map
    inv_map = inp[4]

    print('start generation (ns).', time.time_ns(), 'seed', seed)
    gen_head, gen_tail = naive_sampling(head_index=headss,
                                        tail_index=tailss,
                                        seed=seed)
    print('sample generated (ns).', time.time_ns(), 'seed', seed)
    with open(f"{folder}_Null{seed}_beta10_beta20.tsv", 'w') as out_f:
        if len(inv_map) > 0:
            for idx, h in enumerate(gen_head):
                hs = ','.join([str(inv_map[v]) for v in h])
                ts = ','.join([str(inv_map[v]) for v in gen_tail[idx]])
                out_f.write(f'{hs}\t{ts}\n')
        else:
            for idx, h in enumerate(gen_head):
                hs = ','.join([str(v) for v in h])
                ts = ','.join([str(v) for v in gen_tail[idx]])
                out_f.write(f'{hs}\t{ts}\n')
