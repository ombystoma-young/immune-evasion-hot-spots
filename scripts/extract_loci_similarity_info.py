import os
import argparse
from collections import defaultdict

import pandas as pd
import networkx as nx
from numpy import nan


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_find_novel_adgs" pipeline. '
                                                 'Master script which aggregates information about loci similarity, '
                                                 'extracts ids of desired communities in graph')
    parser.add_argument('-c', '--communities', default=None, type=str, nargs='?',
                        help='path to ClusterONE output')
    parser.add_argument('-g', '--gff', default=None, type=str, nargs='?',
                        help='path to input gff file for target loci')
    parser.add_argument('-e', '--edgelist', default=None, type=str, nargs='?',
                        help='path to edgelist tsv file with information about links between loci (means similarity)')
    parser.add_argument('-m', '--markers', default=None, type=str, nargs='?',
                        help='comma based list of target clans, e.g. "10_clan,106_pair"')
    parser.add_argument('-o', '--outdir', default='.', type=str, nargs='?',
                        help='output directory')
    return parser.parse_args()


# read gff
def _split_attributes(attrs: str):
    pairs = attrs.split(';')
    kvs = [pair.split('=', maxsplit=1) for pair in pairs]
    for i, kv in enumerate(kvs):
        if len(kv) == 1:
            kvs[i - 1][-1] += f'{kv[0]}'
            kvs.remove(kv)
    return {f'ATTRIBUTE_{k}': v for k, v in kvs}


def unites_attributes_col(row: pd.DataFrame, attribute_cols: pd.Index) -> str:
    attributes = []
    for attribute_name in attribute_cols:
        if row[attribute_name] is not nan:
            suffix = '_'.join(attribute_name.split('_')[1:])
            attributes.append(f'{suffix}={row[attribute_name]}')

    return ';'.join(attributes)


def read_gff(gff_path):
    """
    reads gff file into pd dataframe
    :param gff_path: (str) path to gff file
    :return: pd.DataFrame from gff file
    """
    colnames = ['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    entries = []
    with open(gff_path, 'rt') as gff:
        for line in gff:
            entry = line.strip().split('\t')
            entries.append(entry)

    gff_data = pd.DataFrame(entries, columns=colnames)
    gff_data['attribute_dict'] = gff_data['attributes'].apply(_split_attributes)

    norm_attribute = pd.json_normalize(gff_data.attribute_dict)
    gff_data = pd.concat([gff_data, norm_attribute], axis=1)
    # remove temp columns:
    gff_data = gff_data.drop(columns=['attribute_dict', 'attributes'])
    return gff_data


# read info about communities
def get_jaccard_distance(set_1: set, set_2: set) -> float:
    dist = len(set_1 & set_2) / len(set_1 | set_2)
    return dist


# read info from graph
def read_clu1_res(in_path) -> tuple:
    clu2clans = defaultdict(list)
    with open(in_path, 'rt') as in_file:
        for i, line in enumerate(in_file):
            clan = line.strip().split('\t')
            name = f'{i + 1}_community'
            for cluster in clan:
                clu2clans[cluster].append(name)
    return clu2clans, i + 1


def read_edgelist(in_path: str) -> nx.classes.graph.Graph:
    """
    reads the edgefile into the nx object (Graph)
    """
    g = nx.read_edgelist(in_path, delimiter='\t', create_using=nx.Graph, data=(("weight", float),))
    return g


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


def fix_overlaps(clu2clans: dict) -> dict:
    clans = {}
    for cluster in clu2clans.keys():
        if len(clu2clans[cluster]) > 1:
            clans[cluster] = '; '.join([clan_name for clan_name in clu2clans[cluster]])
        else:
            clans[cluster] = clu2clans[cluster][0]
    return clans


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
                    if nei.endswith('community'):
                        structures[node].append(f'outlier: {nei}')
                        break
                    else:
                        structures[node].append('outlier')
            elif all([nei.endswith('community') for nei in outlier_neighbors]):
                if len(outlier_neighbors) == 1:
                    structures[node].append(f'outlier: {outlier_neighbors[0]}')
                else:
                    structures[node].append(f'bridge: {", ".join(set(outlier_neighbors))}')
    return structures


def mark_lonely_clusters(clusters_mono: set, already_known_structures: dict, max_num_prev: int) -> dict:
    """
    mark clusters that does not form any structures (do not form a clans or even pairs), stand alone in the graph
    """
    counter = max_num_prev
    structures = already_known_structures.copy()
    for cluster in clusters_mono:
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
    communities_path = parse_args().communities
    gff_path = parse_args().gff
    edgelist_path = parse_args().edgelist
    markers_str = parse_args().markers
    outdir = parse_args().outdir

    markers = markers_str.split(',')

    # read gff
    gff_df = read_gff(gff_path)
    loci_df = gff_df[['seq_id', 'ATTRIBUTE_clan']]
    loci_df = loci_df.groupby('seq_id').apply(lambda x: frozenset(x['ATTRIBUTE_clan']))
    loci_df = pd.DataFrame(loci_df, columns=['clans'])

    # read result of cl1
    graph = read_edgelist(edgelist_path)
    structures, max_struct = read_clu1_res(communities_path)
    structures, max_struct = find_all_pairs(graph, structures, max_num_prev=max_struct)

    # append singletons
    nodes_in_graph = set(graph.nodes)
    edges_df = pd.DataFrame(graph.edges)
    edges_df['Jaccard_dist'] = edges_df.apply(
        lambda x: get_jaccard_distance(loci_df.loc[x[0]].clans, loci_df.loc[x[1]].clans), axis=1)
    edges_df = pd.concat([edges_df, pd.DataFrame({0: list(set(loci_df.index) - nodes_in_graph)})], ignore_index=True)
    print(f'Number of singletons in loci network: {len(set(loci_df.index) - nodes_in_graph)}')

    # save edgelist with singletons added
    edges_df.to_csv(os.path.join(outdir, 'edge_list_loci_updated.tsv'), index=False,
                    header=False, sep='\t')

    # calculate loci sizes, which will be represented by size of point in network
    loci_df['loc_size'] = loci_df['clans'].apply(lambda x: len(x))
    # save it
    loci_df.to_csv(os.path.join(outdir, 'loci_sizes.tsv'), sep='\t')

    # get info about outliers in graph
    structures = mark_outliers(graph, structures)
    structures = mark_lonely_clusters(set(loci_df.index) - nodes_in_graph, structures, max_struct)

    # save communities into a separate table, will be represented as colors
    structures_wide = fix_overlaps(structures)
    write_clusters(os.path.join(outdir, 'parsed_communities_info.tsv'), structures_wide, fmt='wide')

    # save information about presence of ADGs in loci, will be shown by node shape in the network
    gff_df[['seq_id', 'ATTRIBUTE_clan']][gff_df.ATTRIBUTE_clan.isin(markers)].to_csv(
        os.path.join(outdir, 'has_adgs.tsv'),
        index=False, sep='\t')

    structures_direct = defaultdict(set)
    for seq_id, values in structures.items():
        for structure in values:
            structures_direct[structure].add(seq_id)

    # find communities with desired markers
    has_markers = gff_df[['seq_id', 'ATTRIBUTE_clan']][gff_df.ATTRIBUTE_clan.isin(markers)][
        'seq_id'].to_list()
    structures_with_known = set()
    loci_within_comm_with_known = set()
    for str_name, structure in structures_direct.items():
        subset_loci = loci_df[loci_df.index.isin(structure)]
        subset_loci_markers = subset_loci[subset_loci.index.isin(has_markers)]
        if len(subset_loci_markers) != 0:
            structures_with_known.add(str_name)
            loci_within_comm_with_known = loci_within_comm_with_known | structure
            print(f'{str_name} has known')

    with open(os.path.join(outdir, 'communities_with_known_adgs_htgs.id'), 'wt') as f:
        f.write('\n'.join(loci_within_comm_with_known))
