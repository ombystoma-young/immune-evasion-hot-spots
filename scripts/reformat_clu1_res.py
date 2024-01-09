import argparse
from collections import defaultdict

import networkx as nx


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_process_clustering_results" pipeline. '
                                                 'Parse CL1 results into tsv file')
    parser.add_argument('-i', '--inputtable', default=None,
                        type=str, nargs='?',
                        help='path to input table file')
    parser.add_argument('-e', '--edgelist', default=None,
                        type=str, nargs='?',
                        help='path to input edgelist')
    parser.add_argument('-d', '--descrfile', default=None,
                        type=str, nargs='?',
                        help='path to long table with all possible clusters')
    parser.add_argument('-o', '--outputlong', default=None,
                        type=str, nargs='?',
                        help='path to output table file in long format (for sankey)')
    parser.add_argument('-w', '--outputwide', default=None,
                        type=str, nargs='?',
                        help='path to output table file in wide format (for all other things)')
    return parser.parse_args()


def read_clu1_res(in_path) -> tuple:
    clu2clans = defaultdict(list)
    with open(in_path, 'rt') as in_file:
        for i, line in enumerate(in_file):
            clan = line.strip().split('\t')
            name = f'{i+1}_clan'
            for cluster in clan:
                clu2clans[cluster].append(name)
    return clu2clans, i+1


def read_clu_descr(in_path: str) -> set:
    """
    create a set of clusters from long file
    """
    clusters = set()
    with open(in_path, 'rt') as in_file:
        for line in in_file:
            if not line.startswith('clu'):
                num = line.strip().split('\t')[0]
                clusters.add(f'clu_{num}')
    return clusters


def read_edgelist(in_path: str) -> nx.classes.graph.Graph:
    """
    reads the edgefile into the nx object (Graph)
    """
    g = nx.read_edgelist(in_path, delimiter='\t', create_using=nx.Graph, data=(("weight", float),))
    return g


def fix_overlaps(clu2clans: dict) -> dict:
    clans = {}
    for cluster in clu2clans.keys():
        if len(clu2clans[cluster]) > 1:
            clans[cluster] = '; '.join([clan_name for clan_name in clu2clans[cluster]])
        else:
            clans[cluster] = clu2clans[cluster][0]
    return clans


def find_all_pairs(g, already_known_structures: dict, max_num_prev: int) -> tuple:
    """
    find all pairs of nodes in graph (aka components of size 2)
    """
    structures = already_known_structures.copy()
    counter = max_num_prev

    for component in nx.connected_components(g):
        if len(component) == 2:
            counter += 1
            for node in component:
                name = f'{counter}_pair'
                structures[node].append(name)
    return structures, counter


def mark_outliers(g, already_known_structures: dict) -> dict:
    """
    marks outliers and bridges in graph which are not in communities,
    based ou their neighborhood
    """
    structures = already_known_structures.copy()
    for node in g.nodes():
        if node not in structures.keys():
            outlier_neighbors = []
            for neighbor in g.neighbors(node):
                if neighbor in structures.keys():
                    outlier_neighbors.append(structures[neighbor][0])
                else:
                    outlier_neighbors.append('outlier')

            if 'outlier' in outlier_neighbors:
                for nei in outlier_neighbors:
                    if nei.endswith('clan'):
                        structures[node].append(f'outlier: {nei}')
                        break
                    else:
                        structures[node].append('outlier')
            elif all([nei.endswith('clan') for nei in outlier_neighbors]):
                if len(outlier_neighbors) == 1:
                    structures[node].append(f'outlier: {outlier_neighbors[0]}')
                else:
                    structures[node].append(f'bridge: {", ".join(outlier_neighbors)}')
    return structures


def mark_lonely_clusters(clusters_descr_file: str, already_known_structures: dict, max_num_prev: int) -> dict:
    """
    mark clusters that does not form any structures (do not form a clans or even pairs), stand alone in the graph
    """
    counter = max_num_prev
    structures = already_known_structures.copy()
    clusters = read_clu_descr(clusters_descr_file)
    for cluster in clusters:
        if cluster not in structures.keys():
            counter += 1
            structures[cluster].append(f'{counter}_mono')
    return structures


def write_clusters(out_path: str, clu2clans: dict, fmt: str) -> None:
    if fmt == 'wide':
        with open(out_path, 'wt') as out_file:
            for pair in clu2clans.items():
                out_file.write('\t'.join(pair))
                out_file.write('\n')
    else:
        with open(out_path, 'wt') as out_file:
            for clu, info in clu2clans.items():
                for item in info:
                    out_file.write('\t'.join([clu, item]))
                    out_file.write('\n')


if __name__ == '__main__':
    # parse args
    in_path = parse_args().inputtable
    edgelist_file = parse_args().edgelist
    clusters_descr_file = parse_args().descrfile
    out_long_path = parse_args().outputlong
    out_wide_path = parse_args().outputwide

    # define a graph
    graph = read_edgelist(edgelist_file)

    # get info about clans
    structures, max_struct = read_clu1_res(in_path)

    # get info about pairs in graph
    structures, max_struct = find_all_pairs(graph, structures, max_num_prev=max_struct)

    # get info about outliers in graph
    structures = mark_outliers(graph, structures)

    # get info about lonely clusters (does not for any structure)
    structures = mark_lonely_clusters(clusters_descr_file, structures, max_struct)
    write_clusters(out_long_path, structures, fmt='long')
    # create wide table
    structures_wide = fix_overlaps(structures)
    write_clusters(out_wide_path, structures_wide, fmt='wide')
