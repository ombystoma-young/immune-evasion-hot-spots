import os
import datetime
import argparse
import multiprocessing

import pandas as pd

from math import log10 as lg
from scipy.stats import hypergeom
from numpy.random import permutation, choice


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_find_novel_adgs" pipeline. '
                                                 'Perform vcontact-like calculation of scores for proteins '
                                                 'encoded in permuted (!) loci to find permuted '
                                                 'loci similarity.')
    parser.add_argument('-i', '--input', default=None, type=str, nargs='?',
                        help='path to input gff file')
    parser.add_argument('-n', '--numperm', default=None, type=int, nargs='?',
                        help='Number of permutations to perform')
    parser.add_argument('-t', '--threads', default=None, type=int, nargs='?',
                        help='Number of threads')
    parser.add_argument('-s', '--samplesize', default=None, type=int, nargs='?',
                        help='Size of permuted sample')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output directory, which will contains resulting pickle permutation tables')
    return parser.parse_args()


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


def calculate_total_clan_number(initial_gff_df):
    return len(initial_gff_df['ATTRIBUTE_clan'].unique())


def get_num_of_comparisons(initial_gff_df):
    n = len(initial_gff_df['seq_id'].unique())
    t = n * (n - 1)
    return t


def get_intersect_len(set1, set2):
    return len(set1 & set2)


def get_prob(c, n, k, n_tot):
    """
    Let L_i and L_j be two comparing loci and,
    n := min( |L_i|, |L_j| ) - minimal cardinality of two sets (# clans in smaller set)
    k := max( |L_i|, |L_j| ) - maximal cardinality of two sets (# clans in bigger set)
    c := |L_i intersect L_j |- cardinality of intersection of two sets (# of common clans)
    n_tot := total # of clans
    then we can calculate the probability of obtaining obtaining the result at least 
    as an extreme as the one that was actually observed, given that H0 is true
    (H0: two loci has common clans just by chance).
    The distribution is hypergeometric (X - number of successes c in n independent trials,
    given that there is k successes in population of size n_tot).
    In other words, X - number of clans in smaller set (c), which are also belong to the bigger set,
    given the size of the smaller set (n), size of the bigger set (k), and total number of clans (n_tot)
    """
    return hypergeom.sf(k=c, N=n, n=k, M=n_tot)


def get_score(prob,  mtcc):
    """
    prob := probability of event given H0 (see `get_prob` function),
    mtcc := multiple testing correction coefficient (for Bonferroni correction, 
    controls FWER, see `get_num_of_comparisons` function)
    """
    if prob != 0:
        score = - lg(prob * mtcc)
    else:
        score = float('+inf')
    return score


def _calculate_score(set1, set2, n_total, t):
    # calculate minimal and maximal cavitivity
    n = min((len(set1), len(set2)))
    k = max((len(set1), len(set2)))
    # calculate cavitivity of inte
    c = get_intersect_len(set1, set2)
    pval = get_prob(c=c, n=n, k=k, n_tot=n_total)
    score = get_score(prob=pval, mtcc=t)
    return score


def get_similarity_scores(loci_df):
    n_total = calculate_total_clan_number(loci_df)
    t = get_num_of_comparisons(loci_df)
    # get df of clan sets belonging to each locus
    df = pd.DataFrame(loci_df.groupby('seq_id').apply(lambda x: frozenset(x['ATTRIBUTE_clan'])), 
                      columns=['clans'])
    # calculate similarity score 
    cross_df = pd.DataFrame(df.clans.apply(lambda x: df.clans.apply(lambda y: _calculate_score(x, y, n_total, t))))
    return cross_df


def write_perm_replica(x):
    df, i = x[0].copy(), x[1]
    sample = choice(gff_df.seq_id.unique(), sample_size)
    df = df.query('seq_id in @sample')
    df.loc[:, 'seq_id'] = permutation(df.seq_id.to_list())
    scores_df = get_similarity_scores(df)
    scores_df_melted = scores_df.melt(ignore_index=False)
    scores_df_melted = scores_df_melted[scores_df_melted.seq_id != scores_df_melted.index]
    max_score = scores_df_melted.value.max()
    scores_df_melted.to_pickle(os.path.join(out_dir, f'{i}.pickle'))
    return max_score


if __name__ == '__main__':
    data_path = parse_args().input
    sample_size = parse_args().samplesize
    n_threads = parse_args().threads
    n_permutations = parse_args().numperm
    out_dir = parse_args().output

    os.makedirs(out_dir, exist_ok=False)

    gff_df = read_gff(data_path)
    loci_df = gff_df[['seq_id', 'ATTRIBUTE_clan']]
    args = [(loci_df, i) for i in range(n_permutations)]
    start = datetime.datetime.now()

    with multiprocessing.Pool(processes=n_threads) as pool:
        results = pool.map(write_perm_replica, args)

    n_seconds = datetime.datetime.now() - start
    print(f'{datetime.datetime.now() - start} hours')
