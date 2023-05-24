import os
import re
import json

from collections import OrderedDict

import numpy as np
import pandas as pd


def read_json(json_like: str) -> dict:
    """
    read json string, returns it as a (dict) object
    """
    # based on file structure:
    json_like = json_like.replace(" '", ' "')
    json_like = json_like.replace("{'", '{"')
    json_like = json_like.replace("':", '":')

    js_d = json.loads(json_like)
    return js_d


def define_clu_main_product(products: dict) -> str:
    """

    :param functions: (dict) contains of all functions of cluster
    :return: main_function (str) name of cluster's main functional product
    """
    main_function = None
    max_count = 0
    for product, count in products.items():
        if product != 'hypothetical protein' and count > max_count:
            main_function = product
            max_count = count
    if not main_function and ('hypothetical protein' in products.keys()):
        main_function = 'hypothetical protein'
    return main_function


def extract_cluster_functions_wide(wide_path: str) -> dict:
    """
    extract information about main function of each cluster (dict) {cluster_parent: (number of cluster, function)}
    :param wide_path:
    :return:
    """
    cluster_to_function = {}
    with open(wide_path, 'r') as tsv_file:
        for ind, line in enumerate(tsv_file):
            cluster = line.strip().split('\t')
            clu_parent = cluster[0]
            clu_functs = cluster[2]
            clu_funct_js = read_json(clu_functs)
            clu_main_func = define_clu_main_product(clu_funct_js)
            cluster_to_function[clu_parent] = (ind, clu_main_func)
    return cluster_to_function


def extract_cluster_children_long(long_path: str) -> dict:
    """
    extract information about children in clusters, returns dict {child: parent}
    :param long_path: path to long file with results of upstream gene products clusters
    :return: (dict) {child: parent} in cluster
    """
    childs_to_parent = {}

    with open(long_path, 'r') as tsv_file:
        for line in tsv_file:
            cluster = line.strip().split('\t')
            clu_parent = cluster[0]
            clu_child = cluster[1]
            childs_to_parent[clu_child] = clu_parent

    return childs_to_parent


def read_upstream_gff(gff_path: str) -> dict:
    upstreams = OrderedDict()
    with open(gff_path, 'r') as gff_file:
        for line in gff_file:
            orf_info = line.strip().split('\t')
            chrom = orf_info[0]
            if chrom not in upstreams.keys():
                upstreams[chrom] = [orf_info]
            else:
                upstreams[chrom].append(orf_info)
    return upstreams


def add_cluster_information(upstreams: dict, clu_children: dict, clu_products: dict) -> dict:
    """
    add information about cluster for each cds in gff dict object:
    {genome: [[protein_1 gff], [protein_2 gff]]} -> add_cluster_information ->
    -> {genome: [[protein_1 gff, (num_1, function_1)],
                 [protein_2 gff, (num_2, function_2)]]
        }
    """
    upstreams_mod = upstreams.copy()

    for genome, cdss in upstreams.items():
        for ind, cds in enumerate(cdss):
            attribute = cds[8]
            attribute_list = attribute.split(';')
            cds_id = attribute_list[0].split('=')[1]
            parent = clu_children[cds_id]
            function = clu_products[parent]
            upstreams_mod[genome][ind].append(function)
    return upstreams_mod


def _find_numclu_one(clu_products: dict, keyword: str) -> set:
    """
    :param clu_products (tuple): (num_clu, main_product)
    :param keyword (str): product name, needed to find and exclude
    :return: (set) of cluster numbers to exclude
    """
    tabu_clu_num = set()
    for num, func in clu_products.values():
        if re.findall(keyword, func, re.IGNORECASE):
            tabu_clu_num.add(num)
    return tabu_clu_num


def find_numclu_to_filter(clu_products: dict, keywords: set) -> set:
    """
    uses `_find_numclu_one` for set of keywords
    :param clu_products: clu_products (tuple): (num_clu, main_product)
    :param keywords (set): of product names, needed to find and exclude
    :return: (set) of cluster numbers to exclude
    """
    tabu_clu_nums = set()
    for keyword in keywords:
        tabu_clu_nums = tabu_clu_nums | _find_numclu_one(clu_products, keyword)
    return tabu_clu_nums


def filter_with_conditions(upstream_mod):

    if sum(upstream_mod['strand'] == '-') == len(upstream_mod):  # potential error place
        select_indices = list(np.where(upstream_mod['cluster_in_tabu'])[0])
        target_index = select_indices[0]
    elif sum(upstream_mod['strand'] == '+') == len(upstream_mod):
        select_indices = list(np.where(upstream_mod['cluster_in_tabu'])[0])
        target_index = select_indices[-1]
    else:
        chormosome = upstream_mod.chromosome.unique()[0]
        raise ValueError(f'Non-unique strand for {chormosome}')

    start = upstream_mod.iloc[target_index]['start']
    end = upstream_mod.iloc[target_index]['end']

    if sum(upstream_mod['strand'] == '-') == len(upstream_mod):
        condition = (upstream_mod['end'] <= end) | (upstream_mod['start'] - 5000 >= start)
    else:
        condition = (upstream_mod['end'] >= end) | (upstream_mod['start'] + 5000 <= start)

    upstream_mod = upstream_mod[condition]
    return upstream_mod


def filter_within_one_genome(upstream: list, tabu_clusters: set) -> list:
    """
    finds tabu genes in one genome and returns the resulting list without it's upstream region
    :param upstream (list): part of upstream.gff for one region:
    {genome (list): [[protein_1 gff, (num_1, function_1)],
                 [protein_2 gff, (num_2, function_2)]]
        }
    :param tabu_clu_nums (set): num of tabu
    :return: upstream_mod (list) with removed tabu genes and their upstream
    """
    colnames_ = ['chromosome', 'annotator', 'type', 'start', 'end',
                 'something', 'strand', 'somthing_else', 'attribute', 'cluster']

    upstream_df = pd.DataFrame(upstream, columns=colnames_)  # list to dataframe
    # numbers: str to int
    upstream_df['start'] = list(map(lambda x: int(x), upstream_df['start']))
    upstream_df['end'] = list(map(lambda x: int(x), upstream_df['end']))
    # column: cluster number:
    upstream_df['cluster_num'] = list(map(lambda x: x[0], upstream_df['cluster']))
    # column: is tabu cluster?
    upstream_df['cluster_in_tabu'] = list(map(lambda x: x in tabu_clusters, upstream_df['cluster_num']))

    if sum(upstream_df['cluster_in_tabu']) > 0:
        upstream_df = filter_with_conditions(upstream_df)
    upstream_mod = upstream_df.iloc[:, :10].values.tolist()
    return upstream_mod


def make_gff_with_clusters(path_long: str, path_wide: str, path_upstream: str):
    children_to_clu = extract_cluster_children_long(path_long)
    parent_to_func = extract_cluster_functions_wide(path_wide)
    upstreams = read_upstream_gff(path_upstream)
    ups_with_clu_info = add_cluster_information(upstreams=upstreams,
                                                clu_children=children_to_clu,
                                                clu_products=parent_to_func)
    upstreams_mod = []
    for upstream in ups_with_clu_info.values():
        for cds in upstream:
            upstreams_mod.append(cds)
    return upstreams_mod


def write_modified_gff(upstream_gff: list, gff_path: str) -> None:
    with open(gff_path, 'w') as gff_file:
        for line in upstream_gff:
            line[3] = str(line[3])
            line[4] = str(line[4])
            cluster_num=line[9][0]
            cluster_main_prod=line[9][1]
            gff_file.write('\t'.join(line[:9]))
            gff_file.write(f';{cluster_num=};cluster_main_prod={cluster_main_prod}\n')


if __name__ == '__main__':
    res_dir = 'results'
    prefix = 'upstream_proteins_clu'
    upstream_dir = 'upstream_search'
    path_wide_results = os.path.join(res_dir, f'{prefix}_wide_sorted.tsv')
    path_long_results = os.path.join(res_dir, f'{prefix}_long.tsv')
    upstream_file = os.path.join(upstream_dir, 'upstream_fixed.gff')
    upstream_mod_file = os.path.join(res_dir, 'upstreams_with_clusters.gff')

    modified_upstream = make_gff_with_clusters(path_long_results, path_wide_results, upstream_file)
    write_modified_gff(modified_upstream, upstream_mod_file)