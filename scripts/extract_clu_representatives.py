import argparse
import os

from Bio.SeqIO.FastaIO import SimpleFastaParser


def read_msa_truncated(fasta_path: str) -> list:
    """
    read fasta file into list of tuples and extracts the first one
    :param fasta_path:  (str) path to fasta file with sequence
    :return: tuple (seq_id (str), sequence (str)), corresponding to representative sequence
    """
    with open(fasta_path, 'r') as fasta_file:
        records = list(SimpleFastaParser(fasta_file))
    return records[0]


def extract_representatives(clusters: list, dir_path: str,
                            fasta_path: str) -> list:
    seq_ids = {}
    sequences = []

    for cluster in clusters:
        path = os.path.join(dir_path, cluster)
        seq_id, _ = read_msa_truncated(path)
        seq_ids[seq_id] = cluster
    with open(fasta_path, 'r') as fasta_file:
        for record in SimpleFastaParser(fasta_file):
            if record[0] in seq_ids.keys():
                sequences.append((f'{seq_ids[record[0]]} ({record[0]})', record[1]))
    return sequences


def write_results(sequences: list, out_path: str) -> None:
    """
    writes information about sequences into fasta file
    :param sequences: (list) of (tuples): [(seq_name, sequence), (seq_name, sequence), ...]
    :param out_path: (str) path to output file to be written
    """
    with open(out_path, 'w') as out_file:
        for sequence in sequences:
            out_file.write(f'>{sequence[0]}\n')
            out_file.write(sequence[1])
            out_file.write('\n')


def parse_args():
    parser = argparse.ArgumentParser('Extracts first sequence from query alignments')
    parser.add_argument('-d', '--dir', default=None, nargs='?',
                        help='input directory with alignments')
    parser.add_argument('-f', '--fasta', default=None, nargs='?',
                        help='input file with target sequences')
    parser.add_argument('-s', '--search', default=None, nargs='?', type=str,
                        help='comma-separated list of desired clusters')
    parser.add_argument('-o', '--output', default=None, nargs='?')
    return parser.parse_args()


if __name__ == '__main__':
    dir_path = parse_args().dir
    sequences_path = parse_args().fasta
    query = parse_args().search.split(',')
    output_path = parse_args().output

    sequences = extract_representatives(clusters=query, dir_path=dir_path,
                                        fasta_path=sequences_path)
    write_results(sequences=sequences, out_path=output_path)
