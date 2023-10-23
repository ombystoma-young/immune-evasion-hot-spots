import argparse

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_get_proteomes" pipeline. '
                                                 'Filter predicted by phanotate orfs based on score (-t/--threshold) '
                                                 'and translate sequence')
    parser.add_argument('-i', '--input', default=None, type=str, nargs='?',
                        help='path to raw fasta')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output fasta file')
    parser.add_argument('-t', '--threshold', default=None, type=str, nargs='?',
                        help='threshold of score to consider the orf to be real')
    return parser.parse_args()


def filter_and_translate(in_path: str, t: float) -> list:
    """
    Filter predicted by phanotate orfs based on score and translate sequence.
    :param in_path: (str) path to input file
    :param t: (float) threshold such that values more than t (score > t) are
    considered to be false positive prediction and removed
    :return: (list) of protein sequences
    """
    sequences = []
    with open(in_path, 'rt') as fasta:
        for record in SimpleFastaParser(fasta):
            name = record[0].split()
            score = float(name[2].split('=')[1][:-1])
            if score < t:
                seq = Seq(record[1])
                prot_seq = seq.translate()
                sequences.append((record[0], str(prot_seq)))
    return sequences


def write_fasta(sequences: list, out_path: str) -> None:
    """
    writes list of fasta into file
    :param sequences: (list) of (tuples): (seq name, sequence)
    :param out_path: (str) path to output file
    """
    with open(out_path, 'wt') as out_file:
        for seq in sequences:
            out_file.write('>')
            out_file.write('\n'.join(seq))
            out_file.write('\n')


if __name__ == '__main__':
    input_file = parse_args().input
    threshold = float(parse_args().threshold)
    output = parse_args().output

    seqs = filter_and_translate(in_path=input_file, t=threshold)
    write_fasta(sequences=seqs, out_path=output)