import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser


def get_correspondence(tsv_path: str) -> dict:
    """
    this function extracts relatioships between assembly id and nuccore id
    :param tsv_path: (str) path to tsv file with given information
    :return: (dict) with correspondence: {nuccore_id: assemby_id}
    """
    ids = {}
    with open(tsv_path, 'r') as tsv:
        for line in tsv:
            ids_ = line.strip().split('\t')
            nuccore_id = ids_[0]
            assembly_id = ids_[1]
            ids[nuccore_id] = assembly_id
    return ids


def read_alignments(fasta_path: str) -> list:
    """
    read fasta file into list of tuples
    :param fasta_path:  (str) path to fasta file with sequence
    :return: list of tuples (seq_id (str), sequence (str))
    """
    with open(fasta_path, 'r') as fasta_file:
        records = list(SimpleFastaParser(fasta_file))
    return records


def search_representatives(res_table_path: str, clusters: list) -> dict:
    """
    this function extracts information about representatives from the resulting master-table
    :param res_table_path: path to file, called:
                upstream_proteins_clu_wide_seq_sorted_prokka_domains_pi_length_host.tsv
    :param cluster: (list) of clusters, which representatives are needed
    :return: (dict): {representative assembly_id: cluster_num}
    """
    clu_represents = {}
    with open(res_table_path, 'r') as table:
        for line in table:
            clu = line.strip().split('\t')
            clu_num = clu[14]
            if clu_num != 'cluster_num':
                clu_num = int(clu_num)
            if clu_num in clusters:
                assembly_id = "_".join(clu[0].split('_')[:2])
                clu_represents[assembly_id] = clu_num
    return clu_represents


def get_sequence_of_interest(clu_nums: list, fa_path: str, meta_path: str, table_path: str) -> list:
    """
    extracts particular sequences of cluster representatives from file with MSA, given the clu num.
    :param clu_nums: (list) of clusters of interest
    :param fa_path: (str) path to file with MSA in FASTA format
    :param meta_path: (str) path to file with information about correspondence of nuccore id and assembly id
    :param table_path: (str) resulting file from pipeline of gene findings
    :return: (list) of (tuples) with sequences of cluster representatives
    """
    representatives_seqs = []
    nuccore_to_assembly = get_correspondence(tsv_path=meta_path)
    sequences = read_alignments(fasta_path=fa_path)
    representatives = search_representatives(res_table_path=table_path, clusters=clu_nums)
    representatives_set = set(representatives.keys())
    for sequence in sequences:
        nuccore_id = sequence[0]
        if nuccore_id not in nuccore_to_assembly.keys():
            representatives_seqs.append(sequence)
        elif nuccore_to_assembly[nuccore_id] in representatives_set:
            clu = representatives[nuccore_to_assembly[nuccore_id]]
            representatives_seqs.append((clu, sequence[1]))
    return representatives_seqs


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
    parser = argparse.ArgumentParser()
    parser.add_argument('wildcard', default=None, nargs='?')
    parser.add_argument('input_fa', default=None, nargs='?')
    parser.add_argument('input_meta', default=None, nargs='?')
    parser.add_argument('input_table', default=None, nargs='?')
    parser.add_argument('output_file', default=None, nargs='?')
    return parser.parse_args()


if __name__ == '__main__':
    clusters_nums_ocr = [3, 107, 135, 297]
    clusters_nums_samase = [5, 33, 55, 104, 196, 238, 289, 308, 372, 503, 626]
    clusters = {'ocr': clusters_nums_ocr, 'samase': clusters_nums_samase}
    dataset = parse_args().wildcard
    input_fa = parse_args().input_fa
    input_meta = parse_args().input_meta
    input_table = parse_args().input_table
    output_file = parse_args().output_file
    if input_fa is None:  # manual thing for manual align
        dataset = 'samase'
        input_fa = '../antidefence_trees/upsteam_samase.faa'
        input_meta = '../metadata/assembly_nuccore.tsv'
        input_table = '../results/upstream_proteins_clu_wide_seq_sorted_prokka_domains_pi_length_host.tsv'
        output_file = '../antidefence_trees/representatives_samase.faa'
    sequences = get_sequence_of_interest(clu_nums=clusters[dataset], fa_path=input_fa,
                                         meta_path=input_meta, table_path=input_table)
    write_results(sequences=sequences, out_path=output_file)
