from collections import defaultdict
import sys
import numpy as np


def read_as_edge_list(in_f, vmap):
    # read directed hypergraph as a bipartite graph
    # the suffix 'L' is added to the ids of the hypergraph nodes
    # the suffix 'R' is added to the ids of the hyperedges
    edges = set()
    hmap = dict()
    inv_vmap = {v: k for k, v in vmap.items() if k[0] == "L"}
    for idx, line in enumerate(in_f.readlines()):
        if f'R{idx}' not in hmap:
            nidx = len(hmap)
            hmap[f'R{idx}'] = nidx
        lst = line.strip().split('\t')
        try:
            for i in lst[0].split(','):
                i = i.strip()
                if f'L{i}' not in vmap:
                    nidx = len(vmap)
                    vmap[f'L{i}'] = nidx
                inv_vmap[vmap[f'L{i}']] = f'L{i}'
                edges.add((vmap[f'L{i}'], f'R{idx}'))
        except Exception as e:
            print(e)
            print(lst[0], 'empty head')
        try:
            for i in lst[1].split(','):
                i = i.strip()
                if f'L{i}' not in vmap:
                    nidx = len(vmap)
                    vmap[f'L{i}'] = nidx
                inv_vmap[vmap[f'L{i}']] = f'L{i}'
                edges.add((f'R{idx}', vmap[f'L{i}']))
        except Exception as e:
            print(e)
            print(lst[1], 'empty tail')
    for k, v in hmap.items():
        vmap[k] = len(vmap)
        inv_vmap[vmap[k]] = k
    remapped_edges = list()
    for e in edges:
        if str(e[0])[0] == "R":
            remapped_edges.append((vmap[e[0]], e[1]))
        elif str(e[1])[0] == "R":
            remapped_edges.append((e[0], vmap[e[1]]))
        else:
            print(f"Error in edge {e} remapping.")
            sys.exit(0)
    return inv_vmap, remapped_edges


def read_classes(in_f, vmap):
    # read node labels in a dictionary
    # where the key is the outer node id
    # and the value is the label
    classes = dict()
    dist_classes = dict()
    for line in in_f.readlines():
        lst = line.split('\t')
        c = lst[1].strip()
        if c not in dist_classes:
            dist_classes[c] = len(dist_classes)
        classes[vmap[lst[0].strip()]] = dist_classes[c]
    return classes


def read_as_undirected_hyperedge_list(in_f, vmap):
    # read directed hypergraph as undirected hypergraph
    # vmap stores the outer-inner node ids
    edges = list()
    for line in in_f.readlines():
        lst = line.strip().split('\t')
        tmp = set()
        try:
            for i in lst[0].split(','):
                i = i.strip()
                if i not in vmap:
                    vmap[i] = len(vmap)
            tmp.update([vmap[i.strip()] for i in lst[0].split(',')])
        except Exception as e:
            print(e)
            print(lst[0], 'empty head')
        try:
            for i in lst[1].split(','):
                i = i.strip()
                if i not in vmap:
                    vmap[i] = len(vmap)
            tmp.update([vmap[i.strip()] for i in lst[1].split(',')])
        except Exception as e:
            print(e)
            print(lst[1], 'empty tail')
        edges.append(tuple(sorted(list(tmp))))
    return edges


def read_as_directed_hyperedge_list(in_f, vmap):
    # read directed hypergraph as list of heads and tails
    # vmap stores the outer-inner node ids
    heads = list()
    tails = list()
    for line in in_f.readlines():
        lst = line.strip().split('\t')
        head = set()
        tail = set()
        try:
            for i in lst[0].split(','):
                i = i.strip()
                if i not in vmap:
                    vmap[i] = len(vmap)
            head = [vmap[i.strip()] for i in lst[0].split(',')]
        except Exception as e:
            print(e)
            print(lst[0], 'empty head')
        try:
            for i in lst[1].split(','):
                i = i.strip()
                if i not in vmap:
                    vmap[i] = len(vmap)
            tail = [vmap[i.strip()] for i in lst[1].split(',')]
        except Exception as e:
            print(e)
            print(lst[1], 'empty tail')
        heads.append(tuple(sorted(list(head))))
        tails.append(tuple(sorted(list(tail))))
    return heads, tails


def get_extended_edge_list_from_hypergraph(heads, tails, with_int_edges=True):
    # convert a directed hyperedge to a graph
    # there is an edge between each node in a head and each node in a tail
    # if with_int_edges, then there is an edge between each pair of nodes
    # within a head and each pair of nodes within a tail
    # initialize adjacency matrix
    edge_list = defaultdict(int)
    for idx, h in enumerate(heads):
        for v in h:
            # v can reach each node in the head
            if with_int_edges:
                for w in h:
                    if v != w:
                        edge_list[(v, w)] += 1
            # v can reach each node in the tail
            for w in tails[idx]:
                if v != w:
                    edge_list[(v, w)] += 1
    if with_int_edges:
        for idx, h in enumerate(tails):
            for v in h:
                # v can reach each node in the tail
                for w in h:
                    if v != w:
                        edge_list[(v, w)] += 1
    return [(e[0], e[1], c) for e, c in edge_list.items()]


def read_as_bipartite(in_f, vmap, side='head'):
    # read directed hypergraph as a bipartite graph
    # where right nodes correspond to hyperedges
    # and left nodes correspond to nodes in the head
    # (side='head') or tail (side='tail') of the hyperedges
    edges = list()
    counter = 0
    if side == 'head':
        idx = 0
    else:
        idx = 1
    for line in in_f.readlines():
        lst = line.strip().split('\t')
        for w in range(len(lst)):
            for i in lst[w].split(','):
                i = i.strip()
                if i not in vmap:
                    vmap[i] = len(vmap)
        for i in lst[idx].split(','):
            i = i.strip()
            edges.append((vmap[i], counter))
        counter += 1
    return edges


def process_data_nl_contagion(in_f):
    # read an undirected hypergraph and compute several
    # quantities needed to run the non-linear contagion
    # simulations
    group_list_ = []
    n_list = []
    for line in in_f.readlines():
        hedge = list(set([int(v) for v in line.strip().split(',')]))
        group_list_.append(hedge)
        n_list.append(len(hedge))
    # get node adjacency
    node_mapping = dict()  # relabel nodes from 0 to N-1
    adj_dict = dict()
    for g_id, g in enumerate(group_list_):
        for node in g:
            if node in node_mapping:
                adj_dict[node_mapping[node]].append(g_id)
            else:
                node_id = len(node_mapping)
                node_mapping[node] = node_id
                adj_dict[node_id] = [g_id]
    # relabel group list
    new_group_list = [[node_mapping[node] for node in group]
                      for group in group_list_]
    # get membership
    m_list = [len(adj_dict[node]) for node in adj_dict]
    nmax = max(n_list)
    mmax = max(m_list)
    # get edge-list
    edge_list = []
    for i in adj_dict:
        for j in adj_dict[i]:
            edge_list.append((i, j))
    # get distributions
    pn = np.zeros(nmax + 1)
    for n in n_list:
        pn[n] += 1
    pn /= np.sum(pn)
    gm = np.zeros(mmax + 1)
    for m in m_list:
        gm[m] += 1
    gm /= np.sum(gm)
    # get joint and conditional distributions
    Pmn = np.zeros((nmax + 1, mmax + 1))
    for node in adj_dict:
        for g_id in adj_dict[node]:
            Pmn[n_list[g_id], m_list[node]] += 1
    Pmn /= np.sum(Pmn)  # joint
    Pm_n = np.zeros((nmax + 1, mmax + 1))
    Pn_m = np.zeros((mmax + 1, nmax + 1))
    for n in range(nmax + 1):
        norm = np.sum(Pm_n[n, :])
        if norm > 0:
            Pm_n[n, :] /= norm
    for m in range(mmax + 1):
        norm = np.sum(Pn_m[m, :])
        if norm > 0:
            Pn_m[m, :] /= norm
    results = dict()
    results['group_list'] = new_group_list
    results['adj_dict'] = adj_dict
    results['edge_list'] = edge_list
    results['m_list'] = m_list
    results['n_list'] = n_list
    results['pn'] = pn
    results['gm'] = gm
    results['Pmn'] = Pmn
    results['Pm_n'] = Pm_n
    results['Pn_m'] = Pn_m
    results['nmax'] = nmax
    results['mmax'] = mmax
    return results


def read_and_compute_biot(in_f):
    # read a directed hypergraph and compute
    # its left and right BIOT
    hedges_degs = dict()
    nodes_degs = dict()
    left_edges = []
    right_edges = []
    counter = 0
    for line in in_f.readlines():
        lst = line.strip().split('\t')
        head = [int(v) for v in lst[0].split(',')]
        tail = [int(v) for v in lst[1].split(',')]
        hedges_degs[counter] = (len(head), len(tail))
        for v in head:
            if v not in nodes_degs:
                nodes_degs[v] = [0, 0]
            nodes_degs[v][1] += 1
            left_edges.append((v, counter))
        for v in tail:
            if v not in nodes_degs:
                nodes_degs[v] = [0, 0]
            nodes_degs[v][0] += 1
            right_edges.append((counter, v))
        counter += 1
    left_biot = defaultdict(int)
    for edge in left_edges:
        left_biot[(tuple(nodes_degs[edge[0]]),
                   tuple(hedges_degs[edge[1]]))] += 1
    right_biot = defaultdict(int)
    for edge in right_edges:
        right_biot[(tuple(hedges_degs[edge[0]]),
                    tuple(nodes_degs[edge[1]]))] += 1
    return left_biot, right_biot


def read_and_compute_marginals(in_f):
    # read a directed hypergraph and
    # compute the left and right
    # in- and out-degree sequences
    hedges_degs = dict()
    nodes_degs = dict()
    counter = 0
    for line in in_f.readlines():
        lst = line.strip().split('\t')
        head = [int(v) for v in lst[0].split(',')]
        tail = [int(v) for v in lst[1].split(',')]
        hedges_degs[counter] = (len(head), len(tail))
        for v in head:
            if v not in nodes_degs:
                nodes_degs[v] = [0, 0]
            nodes_degs[v][1] += 1
        for v in tail:
            if v not in nodes_degs:
                nodes_degs[v] = [0, 0]
            nodes_degs[v][0] += 1
        counter += 1
    left_in = sorted([x[0] for x in nodes_degs.values()])
    right_in = sorted([x[0] for x in hedges_degs.values()])
    left_out = sorted([x[1] for x in nodes_degs.values()])
    right_out = sorted([x[1] for x in hedges_degs.values()])
    return left_in, right_in, left_out, right_out
