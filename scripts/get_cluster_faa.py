import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('in_faa', default=None, nargs='?')
    parser.add_argument('in_list', default=None, nargs='?')
    parser.add_argument('out_faa', default=None, nargs='?')
    return parser.parse_args()


def read_tsv(tsv_path: str) -> dict:
    protid2chrom = {}
    with open(tsv_path, 'r') as in_tsv:
        for line in in_tsv:
            entry = line.strip().split('\t')
            chrom = entry[0]
            protein_id = entry[1]
            protid2chrom[protein_id] = chrom
    return protid2chrom


def filter_fasta(protid2chrom: dict, in_path: str) -> dict:
    """
    reads fasta file into dict {seq_id: sequence} with filtering by dict
    :param in_path: (str) path to fasta file to read
    :return: (dict) with correspondence between seq_id and sequence
    """
    sequences = {}
    with open(in_path, 'rt') as fasta:
        for record in SimpleFastaParser(fasta):
            if record[0].split(' ')[0] in protid2chrom.keys():
                sequences[protid2chrom[record[0].split(' ')[0]]] = record[1]
    return sequences


def write_fasta(sequences: dict, out_path: str, delimiter='\n') -> None:
    """
    writes sequences into file in FASTA format
    :param sequences: (dict): {seq_id: sequence}
    :param out_path: (str): path to output file
    :param delimiter: (str): rows delimiter (default: '\n')
    :return: None
    """
    with open(out_path, 'wt') as fasta_out:
        for seqid, sequence in sequences.items():
            fasta_out.write('>')
            fasta_out.write(delimiter.join((seqid, sequence)))
            fasta_out.write(delimiter)


if __name__ == '__main__':
    in_faa = parse_args().in_faa
    in_tsv = parse_args().in_list
    out_faa_ = parse_args().out_faa

    protid2chrom = read_tsv(in_tsv)
    seqs = filter_fasta(protid2chrom, in_faa)
    write_fasta(seqs, out_faa_)
