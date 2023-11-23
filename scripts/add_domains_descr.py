import pprint
import argparse

import pandas as pd
from numpy import nan

COL_NAMES = {'phrog': ['phrog', 'target', 'seqid', 'alnlen',
                       'n_mismatch', 'n_gapopen', 'q_dom_start',
                       'q_dom_end', 't_dom_start', 't_dom_end',
                       'evalue', 'bitscore'],
             'dbapis': ['target_name', 'target_accession',
                        'query_name', 'query_accession',
                        'full_seq_E_value', 'full_seq_score', 'full_seq_bias',
                        'best_one_domain_E_value', 'best_one_domain_score',
                        'best_one_domain_bias',
                        'dom_num_est_exp', 'dom_num_est_reg',
                        'dom_num_est_clu', 'dom_num_est_ov',
                        'dom_num_est_env', 'dom_num_est_dom',
                        'dom_num_est_rep', 'dom_num_est_inc',
                        'description_of_target']
             }
TARGET_MERGE_COLS = {'phrog': 'phrog',
                     'dbapis': 'query_name'}
TARGET_GROUP_COLS = {'phrog': 'target',
                     'dbapis': 'target_name'}
TARGET_SELECT_COLS = {'phrog': ['target', 'phrog', 'annot', 'category'],
                      'dbapis': ['target_name', 'query_name', 'Clan ID', 'APIS genes', 'Defense systems']}


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_process_clustering_results" pipeline. '
                                                 'Read information about clustering by MMSeqs2')
    parser.add_argument('-i', '--input', default=None, type=str, nargs='?',
                        help='path to input table file')
    parser.add_argument('-d', '--description', default=None, type=str, nargs='?',
                        help='path to tsv file with description')
    parser.add_argument('-t', '--type', default=None, type=str, nargs='?',
                        help='path of db used (phrog of dbapis)')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output table file')
    return parser.parse_args()


def read_table(in_path: str, db_type: str) -> pd.DataFrame:
    if db_type == 'phrog':
        df = pd.read_csv(in_path, sep='\t', names=COL_NAMES[db_type])
    else:
        df = pd.read_csv(in_path, sep='\t')
    return df


def fix_na_dbapis(x: pd.DataFrame, type_: str, clans: pd.DataFrame):

    if not x[type_] is nan:
        return x[type_]
    elif x['Clan ID'] == 'CLAN031':
        if type_ == 'Defense systems':
            return 'By homology Thoeris'
        else:
            return 'Homologs Tad1&Tad2'
    elif x['Clan ID'] == 'CLAN018':
        if type_ == 'Defense systems':
            return 'By homology CBASS'
        else:
            return 'Homologs vs.4&Acb2'
    else:
        return clans.loc[x['Clan ID'], type_]


def read_descr(in_path: str, db_type: str) -> pd.DataFrame:
    if db_type == 'phrog':
        annotations = pd.read_csv(in_path, sep='\t')
        annotations['phrog'] = annotations['phrog'].apply(lambda x: f'phrog_{x}')
        annotations['annot'].fillna('hypothetical protein', inplace=True)

    elif db_type == 'dbapis':
        annotations = pd.read_csv(in_path, sep='\t')
        clans = annotations.dropna().set_index('Clan ID')
        annotations['APIS genes'] = (annotations
                                     .apply(lambda x: fix_na_dbapis(x, 'APIS genes', clans), axis=1))
        annotations['Defense systems'] = (annotations
                                          .apply(lambda x: fix_na_dbapis(x, 'Defense systems', clans), axis=1))
        annotations = annotations.rename(columns={'APIS family': 'query_name'}, errors='raise')
    else:
        raise IOError(f'Wrong {db_type=}')
    return annotations


def add_annotation(table: pd.DataFrame, annotation: pd.DataFrame, db_type: str) -> pd.DataFrame:
    table_with_ann = table.merge(annotation, on=TARGET_MERGE_COLS[db_type])
    table_with_ann = table_with_ann[TARGET_SELECT_COLS[db_type]]
    if db_type == 'phrog':
        table_with_ann = (table_with_ann.groupby(TARGET_GROUP_COLS[db_type])
                          .agg(lambda x: ','.join([str(i) for i in set(x)])))
    else:
        table_with_ann = (table_with_ann.groupby(TARGET_GROUP_COLS[db_type])
                          .agg(lambda x: ','.join([str(i) for i in list(x)])))
    return table_with_ann.reset_index()


if __name__ == '__main__':
    in_path = parse_args().input
    descr_path = parse_args().description
    db_type = parse_args().type
    out_path = parse_args().output
    df = read_table(in_path, db_type)
    descr_df = read_descr(descr_path, db_type)
    out_df = add_annotation(table=df,
                            annotation=descr_df,
                            db_type=db_type)
    out_df.to_csv(out_path, sep='\t', index=False)
