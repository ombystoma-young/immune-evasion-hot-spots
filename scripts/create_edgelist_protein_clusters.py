import argparse
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_process_clustering_results" pipeline. '
                                                 'Creates an edgefile from hhsearch blasttbale resulting file')
    parser.add_argument('-i', '--input', default=None, type=str, nargs='?',
                        help='path to input table file')
    parser.add_argument('-p', '--prob', default=None, type=float, nargs='?',
                        help='min prob considered as true')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output table file')
    return parser.parse_args()


def create_edgelist(source_df: pd.DataFrame, prob: float) -> pd.DataFrame:
    df = source_df.copy()
    df['edge_pair'] = df.apply(lambda row: frozenset((row['Query'], row['Hit'])), axis=1)
    df = df.drop(columns=['Query', 'Hit'])
    df_max_prob_filt = df.groupby(by=['edge_pair']).max()
    df_max_prob_filt = df_max_prob_filt.query('`Prob` > @prob')
    df_max_prob_filt['weight'] = df_max_prob_filt['Prob'] / 100
    df_max_prob_filt = df_max_prob_filt.reset_index()
    df_max_prob_filt = df_max_prob_filt[df_max_prob_filt['edge_pair'].map(len) != 1]
    df_max_prob_filt['edge_pair_list'] = df_max_prob_filt['edge_pair'].apply(lambda x: list(x))
    df_max_prob_filt['first'] = df_max_prob_filt['edge_pair_list'].apply(lambda x: f'clu_{int(x[0])}')
    df_max_prob_filt['second'] = df_max_prob_filt['edge_pair_list'].apply(lambda x: f'clu_{int(x[1])}')
    df_edgelist = df_max_prob_filt[['first', 'second', 'weight']]
    return df_edgelist


if __name__ == '__main__':
    in_path = parse_args().input
    prob = parse_args().prob
    out_path = parse_args().output
    hh_res_df = pd.read_csv(in_path, sep='\t', names=['No', 'Hit', 'Prob', 'E-value',
                                                      'P-value', 'Score', 'SS', 'Cols',
                                                      'Length', 'Query_start', 'Query_end',	'Template_start', 'Template_end', 'Query'])
    hh_res_df = hh_res_df[['Query', 'Hit', 'Prob']]
    edge_list = create_edgelist(hh_res_df, prob)
    edge_list.to_csv(out_path, sep='\t', index=False, header=False)
