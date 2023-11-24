import argparse
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_process_clustering_results" pipeline. '
                                                 'Creates an edgefile from hhsearch blasttbale resulting file')
    parser.add_argument('-i', '--input', default=None, type=str, nargs='?',
                        help='path to input table file')
    parser.add_argument('-e', '--ethres', default=None, type=float, nargs='?',
                        help='max evalue considered as true')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output table file')
    return parser.parse_args()


def create_edgelist(source_df: pd.DataFrame, evalthres: float) -> pd.DataFrame:
    df = source_df.copy()
    df['edge_pair'] = df.apply(lambda row: frozenset((row['query'], row['target'])), axis=1)
    df = df.drop(columns=['query', 'target'])
    df_min_eval = df.groupby(by=['edge_pair']).min()
    df_min_eval_filt = df_min_eval.query('`eval` < @evalthres')
    df_min_eval_filt['weight'] = - np.log10(df_min_eval_filt['eval'])
    df_min_eval_filt.replace(np.inf, 116, inplace=True)
    df_min_eval_filt = df_min_eval_filt.reset_index()
    df_min_eval_filt = df_min_eval_filt[df_min_eval_filt['edge_pair'].map(len) != 1]
    df_min_eval_filt['edge_pair_list'] = df_min_eval_filt['edge_pair'].apply(lambda x: list(x))
    df_min_eval_filt['first'] = df_min_eval_filt['edge_pair_list'].apply(lambda x: f'clu_{int(x[0])}')
    df_min_eval_filt['second'] = df_min_eval_filt['edge_pair_list'].apply(lambda x: f'clu_{int(x[1])}')
    df_edgelist = df_min_eval_filt[['first', 'second', 'weight']]
    return df_edgelist


if __name__ == '__main__':
    in_path = parse_args().input
    e_value = parse_args().ethres
    out_path = parse_args().output
    hh_res_df = pd.read_csv(in_path, sep='\t', names=['query', 'target', 'match/tLen', 'eval'])
    hh_res_df = hh_res_df[['query', 'target', 'eval']]
    edge_list = create_edgelist(hh_res_df, e_value)
    edge_list.to_csv(out_path, sep='\t', index=False, header=False)
