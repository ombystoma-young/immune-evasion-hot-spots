import os
import numpy as np
import pandas as pd

from math import log10 as lg
from scipy.stats import hypergeom
from tqdm import tqdm

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


# calculate similarity
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


def get_score(prob, mtcc):
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


if __name__ == '__main__':
    max_perm_values = []
    n_permutations = 5000
    data_path = 'data_autographiviridae_refseq/upstreams/early_with_clusters.gff'
    in_dir = 'data_autographiviridae_refseq/loci_similarity/permutations'
    out_file = 'data_autographiviridae_refseq/loci_similarity/edge_list_loci.tsv'
    # calculate scores

    gff_df = read_gff(data_path)
    loci_df = gff_df[['seq_id', 'ATTRIBUTE_clan']]
    scores_df = get_similarity_scores(loci_df)
    scores_df_melted = scores_df.melt(ignore_index=False)
    scores_df_melted = scores_df_melted[scores_df_melted.seq_id != scores_df_melted.index]

    # calculate max value from permutations
    for i in tqdm(range(n_permutations)):
        scores_df_melted_perm = pd.read_pickle(os.path.join(in_dir, f'{i}.pickle'))
        scores_df_melted_perm = scores_df_melted_perm[scores_df_melted_perm.value != float('+inf')]
        max_perm_values.append(scores_df_melted_perm.value.max())
    # write results

    max_perm_values = np.array(max_perm_values)
    scores_df_melted[scores_df_melted.value >= max_perm_values.max()].loc[:, 'seq_id'].to_csv(out_file,
                                                                                              sep='\t', header=False)
