import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_meta_annotate_genomes" pipeline. '
                                                 'Parse information about completeness,'
                                                 ' obtained with CheckV completeness, and '
                                                 'finds contigs with high completeness')
    parser.add_argument('-i', '--input', default=None, type=str, nargs='?',
                        help='path to input tsv file')
    parser.add_argument('-t', '--threshold', default=90.0, type=float, nargs='?',
                        help='threshold to consider the contig to be complete and take it to '
                             'further analysis')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output txt file with ids of contigs')
    return parser.parse_args()


def read_table(in_path: str) -> pd.DataFrame:
    """
    reads CheckV results table
    :param in_path: (path) to input file
    :return:
    """
    df = pd.read_csv(in_path, sep='\t')
    return df


def parse_completeness_one_contig(row: pd.DataFrame, thres: float) -> bool:
    if row['aai_confidence'] == 'low':
        score = row['hmm_completeness_lower']
    else:
        score = row['aai_completeness']
    return score >= thres


def extract_complete_contigs(source_df: pd.DataFrame, thres: float) -> pd.DataFrame:
    df = source_df.copy()
    df['is_complete'] = df.apply(lambda row:
                                 parse_completeness_one_contig(row,
                                                               thres), axis=1)
    return df[df['is_complete']]


if __name__ == '__main__':
    in_path = parse_args().input
    thres = parse_args().threshold
    out_path = parse_args().output
    df = read_table(in_path)
    df = extract_complete_contigs(df, thres)
    df['contig_id'].to_csv(out_path, header=False, index=False)
