import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_filter_early" pipeline. '
                                                 'Finds length of nucleotide sequences and write them to .bed file: '
                                                 'chrm \t 1 \t length')
    parser.add_argument('-f', '--fna', default=None, type=str, nargs='?',
                        help='path to fasta')
    parser.add_argument('-b', '--bed', default=None, type=str, nargs='?',
                        help='path to output bed')
    return parser.parse_args()


def get_chr_length_couples(in_path: str) -> list:
    """
    reads fasta file into dict {seq_id: sequence}
    :param in_path: (str) path to fasta file to read
    :return: (list) of tuples (seq_id, seq_len)
    """
    sequences = []
    with open(in_path, 'rt') as fasta:
        for record in SimpleFastaParser(fasta):
            seq_id = record[0].split()[0]
            seq_len = len(record[1])
            sequences.append((seq_id, seq_len))
    return sequences


def write_bed(couples_list: list, bed_name: str) -> None:
    with open(bed_name, 'wt') as bed_file:
        for couple in couples_list:
            chrm = couple[0]
            length = couple[1]
            line = "\t".join([chrm, str(0), str(length)])
            bed_file.write(line)
            bed_file.write('\n')


if __name__ == '__main__':
    in_fasta = parse_args().fna
    out_bed = parse_args().bed
    couples = get_chr_length_couples(in_fasta)
    write_bed(couples_list=couples, bed_name=out_bed)
