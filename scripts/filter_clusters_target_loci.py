import argparse

import pandas as pd
from unite_protein_info_to_clusters import group_by_clu


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_find_novel_adgs" pipeline. '
                                                 'Filter long table and calculate table for '
                                                 'selected loci only')
    parser.add_argument('--long', default='../data_autographiviridae_meta/clans_info/res_table_long.tsv', type=str, nargs='?',
                        help='path to input table file without grouping by clusters')
    parser.add_argument('--ids', default='../data_autographiviridae_meta/loci_similarity/communities_with_known_adgs_htgs.id', type=str, nargs='?',
                        help='path to input txt file with selected_ids')
    parser.add_argument('--output', default=None, type=str, nargs='?',
                        help='path to output table file')

    return parser.parse_args()


def read_ids(in_path: str) -> set:
    records = set()
    with open(in_path, 'rt') as in_file:
        for line in in_file:
            records.add(line.strip())
    return records


if __name__ == '__main__':
    tab_long = parse_args().long
    ids_file = parse_args().ids
    res_tab = parse_args().output

    ids = read_ids(ids_file)
    agg_prot_df = pd.read_csv(tab_long, sep='\t')
    # filter long table
    agg_prot_df['seq_id'] = agg_prot_df['prot'].str.split('.').str[:-1].str.join('.')
    agg_prot_df = agg_prot_df[agg_prot_df['seq_id'].isin(ids)]
    agg_prot_df = agg_prot_df.fillna('nan')
    # aggregate clusters results
    agg_grouped_prot_df = group_by_clu(agg_prot_df)
    agg_grouped_prot_df['clu_netw'] = agg_grouped_prot_df['clu'].apply(lambda x: f'clu_{x}')
    agg_grouped_prot_df.to_csv(res_tab, sep='\t', index=False)
