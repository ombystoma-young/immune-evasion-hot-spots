import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_filter_early" pipeline. '
                                                 'Extracts target CDSs from the file with all CDSs')
    parser.add_argument('-i', '--input', default=None, type=str, nargs='?',
                        help='path to input protein fasta with all CDSs')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output faa filte of interest')
    parser.add_argument('-g', '--gff', default=None, type=str, nargs='?',
                        help='path to gff with CDSs of interest')
    return parser.parse_args()


def read_gff(path_to_gff):
    cds_ids = []
    with open(path_to_gff, 'r') as gff:
        for line in gff:
            cds = line.strip().split('\t')
            attribute = cds[-1]
            attribute_list = attribute.split(';')
            cds_id = attribute_list[0].split('=')[1]
            cds_ids.append(cds_id)
    return set(cds_ids)


def read_fasta(in_path: str) -> dict:
    """
    reads fasta file into dict {seq_id: sequence}
    :param in_path: (str) path to fasta file to read
    :return: (dir) with correspondence between seq_id and sequence
    """
    sequences = {}
    with open(in_path, 'rt') as fasta:
        for record in SimpleFastaParser(fasta):
            name = record[0].split()[0]
            sequences[name] = record[1]
    return sequences


def find_upstream_faa(gff_path, faa_input_path, faa_output_path):
    seqs = read_fasta(faa_input_path)
    cds_ids = read_gff(gff_path)
    with open(faa_output_path, 'w') as out_faa:
        for cds_id in seqs.keys():
            if cds_id in cds_ids:
                out_faa.write(f'>{cds_id}\n')
                out_faa.write(f'{seqs[cds_id]}\n')


if __name__ == '__main__':
    in_file = parse_args().input
    out_file = parse_args().output
    gff = parse_args().gff

    find_upstream_faa(gff_path=gff,
                      faa_input_path=in_file,
                      faa_output_path=out_file)
