from copy import deepcopy
import heapq
import zipfile
import io
from tqdm import tqdm
import numpy as np
import sys
sys.path.insert(1, '../')
from helpers.parser_for_output_file_names import parser_for_output_file_name as parse_file_name
from helpers.io import read_as_directed_hyperedge_list
from helpers.utils import compute_v_hlist_maps
from helpers.config import met_path


def compute_m_io_degrees(m, vhead_map, vtail_map, hsize_map):
    '''
    INPUT
    ======
    m, int: min size of a hyperedge to be considered
            in the computation of the degrees
    vhead_map, dict: for each vertex v, list of hyperedges with v in the head
    vtail_map, dict: for each vertex v, list of hyperedges with v in the tail
    hsize_map, dict: size of each hyperedge

    OUTPUT
    ======
    io_degrees, dict: in- and out-degree of each vertex,
                      restricted to hyperedges of size >= m
    m_vhead_map, dict: for each vertex v, list of hyperedges of
                       size >= m with v in the head
    m_vtail_map, dict: for each vertex v, list of hyperedges of
                       size >= m with v in the tail
    '''
    io_degrees = dict()
    m_vhead_map = dict()
    m_vtail_map = dict()
    # get all distinct vertices

    vertices = set(vhead_map.keys())
    vertices.update(vtail_map.keys())
    print('vertices', len(vertices))
    for v in vertices:
        ideg = 0
        odeg = 0
        for h in vhead_map.get(v, []):
            if hsize_map.get(h, 0) >= m:
                odeg += 1
                lst = m_vhead_map.get(v, [])
                lst.append(h)
                m_vhead_map[v] = lst
        for h in vtail_map.get(v, []):
            if hsize_map.get(h, 0) >= m:
                ideg += 1
                lst = m_vtail_map.get(v, [])
                lst.append(h)
                m_vtail_map[v] = lst
        io_degrees[v] = (ideg, odeg)
    return io_degrees, m_vhead_map, m_vtail_map


def compute_m_coreness(m, n, degrees,
                       vhead_map,
                       vtail_map,
                       hsize_map,
                       heads, tails,
                       side):
    """
    INPUT
    ======
    m, int: min size of a hyperedge to be considered
            in the computation of the degrees
    n, int: number of vertices in the hypergraph
    degrees, dict: in- and out-degree of each vertex,
                      restricted to hyperedges of size >= m
    vhead_map, dict: for each vertex v, list of hyperedges of
                     size >= m with v in the head
    vtail_map, dict: for each vertex v, list of hyperedges of
                     size >= m with v in the tail
    hsize_map, dict: size of each hyperedge
    heads, list: heads of the hyperedges
    tails, list: tails of the hyperedges
    side, int: if 0 we peel the hypergraph based on the in-degree;
               if 1 we peel based on the out-degree
    OUTPUT
    ======
    core_values, np.array: m-i-coreness (side = 0) or m-o-coreness (side = 1)
                           of each vertex
    """
    coreness = np.zeros(n, dtype=np.int32)
    max_deg = max([i[side] for i in degrees.values()])
    for k in range(max_deg):
        # initialize degrees
        io_degrees = deepcopy(degrees)
        k_core = set()
        # vertices visited
        visited = set()
        # hyperedges removed
        removed = set()
        # initialize heap
        mh = []
        for v, (i, o) in io_degrees.items():
            if side == 0:
                mh.append((i, v))
            elif side == 1:
                mh.append((o, v))
        heapq.heapify(mh)
        while len(mh) > 0 and len(visited) < n:
            # current node of minimum degree
            _, min_node = heapq.heappop(mh)
            while min_node in visited and len(mh) > 0:
                _, min_node = heapq.heappop(mh)
            if min_node in visited and len(mh) == 0:
                break
            visited.add(min_node)
            if io_degrees[min_node][side] >= k:
                k_core.add(min_node)
            else:
                # decrease the size of the hyperedges
                # to simulate the deletion of min_node
                h_changed = set(vhead_map.get(min_node, []))
                h_changed.update(vtail_map.get(min_node, []))
                for h in h_changed:
                    hsize_map[h] -= 1
                # if the hyperedge has now size < m, we decrease the degree
                # of the nodes in such hyperedge to simulate its deletion
                # contains the hyperedges that will be removed
                for h in h_changed:
                    if hsize_map.get(h, 0) < m:
                        if h in removed:
                            continue
                        if side == 1:
                            for v in heads[h]:
                                if v not in visited:
                                    i, o = io_degrees.get(v, (0, 0))
                                    io_degrees[v] = (i, o - 1)
                                    heapq.heappush(mh, (o - 1, v))
                        elif side == 0:
                            for v in tails[h]:
                                if v not in visited:
                                    i, o = io_degrees.get(v, (0, 0))
                                    io_degrees[v] = (i - 1, o)
                                    heapq.heappush(mh, (i - 1, v))
                        removed.add(h)
        for v in k_core:
            coreness[v] = k
        if len(k_core) == 0:
            break
    return coreness


def compute_m_io_corenesses(max_m, n,
                            vhead_map, vtail_map,
                            hsize_map,
                            heads, tails, file_path):
    """
    INPUT
    ======
    max_m, int: max m-coreness to compute
    n, int: number of vertices in the hypergraph
    vhead_map, dict: for each vertex v, list of hyperedges of
                     size >= m with v in the head
    vtail_map, dict: for each vertex v, list of hyperedges of
                     size >= m with v in the tail
    hsize_map, dict: size of each hyperedge
    heads, list: heads of the hyperedges
    tails, list: tails of the hyperedges
    file_path, str: path to store output files
    """
    out_f_i = open(file_path + '_i.tsv', "w")
    out_f_o = open(file_path + '_o.tsv', "w")
    coreness_header = 'VertexId\tm\tm-shellIndex\n'
    out_f_i.write(coreness_header)
    out_f_o.write(coreness_header)
    for m in tqdm(range(2, max_m)):
        print('examining', m)
        iodeg, mvhead, mvtail = compute_m_io_degrees(m,
                                                     vhead_map,
                                                     vtail_map,
                                                     hsize_map)
        hsize_map_2 = deepcopy(hsize_map)
        iodeg_2 = deepcopy(iodeg)
        icoreness = compute_m_coreness(m,
                                       n,
                                       iodeg,
                                       mvhead,
                                       mvtail,
                                       hsize_map,
                                       heads,
                                       tails,
                                       0)
        for vidx, c in enumerate(icoreness):
            out_f_i.write(f'{vidx}\t{m}\t{c}\n')
        print('icoreness computed.')
        ocoreness = compute_m_coreness(m,
                                       n,
                                       iodeg_2,
                                       mvhead,
                                       mvtail,
                                       hsize_map_2,
                                       heads,
                                       tails,
                                       1)
        for vidx, c in enumerate(ocoreness):
            out_f_o.write(f'{vidx}\t{m}\t{c}\n')
        print('ocoreness computed.')
        if np.sum(icoreness) == 0 and np.sum(ocoreness) == 0:
            break
    out_f_i.close()
    out_f_o.close()
    return


def compute_hyper_coreness_from_array(shell_indices):
    '''
    INPUT
    ======
    shell_indices, dict: for each m, the array of m-shell values of each v

    OUTPUT
    ======
    hypercoreness, dict: for each vertex, the hyper-coreness
    '''
    hypercoreness = {v: 0. for v in range(len(shell_indices[0]))}
    for m in range(2, len(shell_indices)):
        max_coreness = max(shell_indices[m])
        if max_coreness == 0:
            max_coreness = 1
        # normalize coreness values by max
        for v in range(len(shell_indices[m])):
            hypercoreness[v] += (shell_indices[m][v] / max_coreness)
    return hypercoreness


def run_coreness(inp):
    '''
    Computes the m-i-shell index and the m-o-shell index of each vertex for
    each hyperedge size m. The m-i-shell index is the in-degree of the vertex
    when it is removed from the hypergraph during the in-degree based peeling.
    The m-o-shell index is the out-degree of the vertex when it is removed
    from the hypergraph during the out-degree based peeling.
    Finally, the i-hyper-coreness (o-hyper-coreness) of a vertex is the sum
    of its m-i-shell (m-o-shell) indices.
    '''
    file_path = inp[0]
    file_name = inp[1]
    vmap = inp[2]
    max_m = inp[3]
    num_vertices = inp[4]
    root = inp[5]
    typ = inp[6]

    if file_path.endswith('.zip'):
        z = zipfile.ZipFile(file_path, "r")
        with z.open(file_name, 'r') as f1:
            map_fields = parse_file_name(file_name, typ)
            fobj = io.TextIOWrapper(f1, encoding='utf-8', newline='')
            sheads, stails = read_as_directed_hyperedge_list(fobj, vmap)
    else:
        with open(file_path) as f1:
            map_fields = dict()
            sheads, stails = read_as_directed_hyperedge_list(f1, vmap)

    v_to_headlist, v_to_taillist, hsize_map = compute_v_hlist_maps(sheads,
                                                                   stails)

    if "randomSeed" in map_fields:
        core_path = f'{met_path}{root}__randomSeed_{map_fields["randomSeed"]}'
        core_path += f'__algorithm_{map_fields["algorithm"]}__coreness'
    else:
        core_path = f'{met_path}{root}_coreness'

    compute_m_io_corenesses(max_m,
                            num_vertices,
                            v_to_headlist,
                            v_to_taillist,
                            hsize_map,
                            sheads,
                            stails,
                            core_path)
