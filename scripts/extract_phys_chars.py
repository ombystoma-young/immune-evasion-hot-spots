import argparse
import pandas as pd

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_process_clustering_results" pipeline. '
                                                 'Calculates physical properties of AA sequences in fasta file')
    parser.add_argument('-i', '--input', default=None, type=str, nargs='?',
                        help='path to input table file')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output table file')
    return parser.parse_args()


def read_sequences(in_path: str) -> pd.DataFrame:
    with open(in_path, 'rt') as fasta_file:
        records = list(SimpleFastaParser(fasta_file))

    prot_seq_df = pd.DataFrame(records, columns=['protein_ID', 'sequence'])
    return prot_seq_df


def calculate_pi(seq: str) -> float:
    pi = IP(seq).pi()
    return round(pi, 2)


def calculate_seq_len(seq: str) -> int:
    return len(seq)


def calculate_phys_char(source_df: pd.DataFrame) -> pd.DataFrame:
    df = source_df.copy()
    df['length'] = df['sequence'].apply(lambda seq: calculate_seq_len(seq))
    df['pi'] = df['sequence'].apply(lambda seq: calculate_pi(seq))
    return df


if __name__ == '__main__':
    in_file = parse_args().input
    out_file = parse_args().output
    seqs_df = read_sequences(in_file)
    seqs_df = calculate_phys_char(seqs_df)
    seqs_df.to_csv(out_file, sep='\t', header=False, index=False)
