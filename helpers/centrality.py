import sys
sys.path.insert(1, '../')
from helpers.io import read_as_directed_hyperedge_list, get_extended_edge_list_from_hypergraph, read_as_edge_list
from helpers.parser_for_output_file_names import parser_for_output_file_name as parse_file_name
from copy import deepcopy
from sknetwork.ranking import PageRank, HITS
from sknetwork.data import parse
import io
import zipfile


def parallel_pagerank(inp):
    '''
    Run PageRank for the sample read from file_name.
    '''
    file_path = inp[0]
    file_name = inp[1]
    vmap = inp[2]
    typ = inp[3]

    z = zipfile.ZipFile(file_path, "r")
    with z.open(file_name, 'r') as f1:
        map_fields = parse_file_name(file_name, typ)
        fobj = io.TextIOWrapper(f1, encoding='utf-8', newline='')
        sheads, stails = read_as_directed_hyperedge_list(fobj, vmap)
    sext_edge_list = get_extended_edge_list_from_hypergraph(sheads,
                                                            stails,
                                                            with_int_edges=False)
    sAdj = parse.from_edge_list(sext_edge_list,
                                reindex=False,
                                matrix_only=True,
                                directed=True,
                                weighted=True)
    pagerank = PageRank(n_iter=5, damping_factor=0.85,
                        solver='piteration', tol=1e-6)

    sscores_pr = pagerank.fit_predict(sAdj)
    # Node Id, Sampler, Sample, PageRank
    rows = []
    for idx, pr in enumerate(sscores_pr):
        rows.append([str(idx),
                     map_fields['algorithm'],
                     str(map_fields['randomSeed']),
                     str(pr)])
    return rows


def map_node_to_type(x):
    if x[0] == 'L':
        return 'Vertex'
    if x[0] == 'R':
        return 'Hyperedge'
    return "None"


def parallel_hits(inp):
    '''
    Run HITS for the sample read from file_name.
    '''
    file_path = inp[0]
    file_name = inp[1]
    vmap = inp[2]
    typ = inp[3]

    vmap2 = deepcopy(vmap)
    z = zipfile.ZipFile(file_path, "r")
    with z.open(file_name, 'r') as f1:
        map_fields = parse_file_name(file_name, typ)
        fobj = io.TextIOWrapper(f1, encoding='utf-8', newline='')
        inv_vmap, sdir_edges = read_as_edge_list(fobj, vmap2)
    # in the case of UnpretentiousNullModel, some nodes in the original
    # hypergraph do not exist in the samples, and thus we treat them as
    # disconnected nodes
    for v, k in vmap.items():
        if v[0] == "L" and k not in inv_vmap:
            inv_vmap[k] = v
    sadj_mat = parse.from_edge_list(sdir_edges,
                                    reindex=False,
                                    matrix_only=True,
                                    directed=True,
                                    bipartite=False,
                                    weighted=False)
    # run HITS
    hits = HITS()
    hits.fit(sadj_mat)
    sscores_h, sscores_a = hits.scores_row_, hits.scores_col_

    # Node Id, Sampler, Sample, Hub Score, Authority Score, Type
    rows = []
    for idx, hs in enumerate(sscores_h):
        rows.append([str(idx),
                     map_fields['algorithm'],
                     str(map_fields['randomSeed']),
                     str(hs),
                     str(sscores_a[idx]),
                     map_node_to_type(inv_vmap.get(idx, "None"))])
    return rows
