import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def compute_statistics(data, id_head, id_tail):
    '''
    Compute characteristics of the directed hypergraph
    in 'data'. Each row in 'data' is a tuple corresponding
    to a directed hyperedge: 'id_head' is the position of the
    head in the tuple, and 'id_tail' is the position of the tail
    in the tuple. 
    We assume that head and tail are comma-separated strings of 
    node ids.
    '''
    vertices = set()
    head_vmap = dict()
    tail_vmap = dict()
    # head, tail, total
    sizes = []
    for row in data:
        head = set(row[id_head].strip().split(','))
        tail = set(row[id_tail].strip().split(','))
        sizes.append([len(head), len(tail), len(head.union(tail))])
        for v in head:
            v = v.strip()
            vertices.add(v)
            head_vmap[v] = head_vmap.get(v, 0) + 1
        for v in tail:
            v = v.strip()
            vertices.add(v)
            tail_vmap[v] = tail_vmap.get(v, 0) + 1
    sizes_df = pd.DataFrame(sizes, columns=['Head',
                                            'Tail',
                                            'Total'], dtype=np.int32)
    print(sizes_df.describe().round(4))
    avg_h = round(1.0 * sum(head_vmap.values()) / len(head_vmap), 4)
    avg_t = round(1.0 * sum(tail_vmap.values()) / len(tail_vmap), 4)
    print(f'Num Vertices: {len(vertices)}')
    print(f'Num HyperEdges: {len(data)}')
    print(f'Avg Out Degree: {avg_h}')
    print(f'Avg In Degree: {avg_t}')
    # plot
    fig, ax = plt.subplots(1, 2, figsize=(20, 6))
    sns.histplot(data=sizes_df, x='Head', kde=True, ax=ax[0])
    sns.histplot(data=sizes_df, x='Tail', kde=True, ax=ax[1])


def compute_v_hlist_maps(heads, tails, size=1):
    '''
    INPUT
    ======
    heads, list: list of hyperedge heads
    tails, list: list of hyperedge tails
    size, int: min size for a hyperedge to be considered

    OUTPUT
    ======
    v_to_hlist, dict: for each vertex, list of hyperedge heads of size >= size
                      including it
    v_to_tlist, dict: for each vertex, list of hyperedge tails of size >= size
                      including it
    h_to_size, dict: for each hyperedge of size >= size, its size
    '''
    v_to_hlist = dict()
    v_to_tlist = dict()
    h_to_size = dict()
    for idx, h in enumerate(heads):
        tmp = set(h)
        tmp.update(tails[idx])
        if len(tmp) < size:
            continue
        h_to_size[idx] = len(tmp)
        for v in h:
            # update v map
            lst = v_to_hlist.get(v, set())
            lst.add(idx)
            v_to_hlist[v] = lst
        for v in tails[idx]:
            # update v map
            lst = v_to_tlist.get(v, set())
            lst.add(idx)
            v_to_tlist[v] = lst
    return v_to_hlist, v_to_tlist, h_to_size
