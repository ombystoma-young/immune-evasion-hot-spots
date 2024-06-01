import os
import json
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_extract_taxonomy_refseq" pipeline. '
                                                 'Finds all taxon ids for all downloaded sequences')
    parser.add_argument('-i', '--input_dir', default=None, type=str, nargs='?',
                        help='path to NCBI Datasets assembly report')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output fasta file')
    return parser.parse_args()


def extract_jsonl(genomes_dir: str) -> dict:
    assembly_to_host = {}
    nuccore_to_host = {}
    file_path = os.path.join(genomes_dir, 'assembly_data_report.jsonl')
    with open(file_path, 'r') as jsonl:
        for line in jsonl:
            entry = json.loads(line.strip())
            organism_name = entry['organism']['organismName']
            host = organism_name.split()[0]
            assembly_id = entry['accession']
            assembly_to_host[assembly_id] = host
    for assembly_id in assembly_to_host.keys():
        seq_rep = os.path.join(genomes_dir, assembly_id, 'sequence_report.jsonl')
        with open(seq_rep, 'rt') as in_file:
            for line in in_file:
                entry = json.loads(line.strip())
                if assembly_id.startswith('GCA'):
                    nuccore = entry['genbankAccession']
                else:
                    nuccore = entry['refseqAccession']
                nuccore_to_host[nuccore] = assembly_to_host[assembly_id]
    return nuccore_to_host


if __name__ == '__main__':
    gen_dir = parse_args().input_dir
    output = parse_args().output

    nuccore_to_host = extract_jsonl(genomes_dir=gen_dir)

    with open(output, 'w') as tsv:
        for id_ in nuccore_to_host.keys():
            tsv.write(f'{id_}\t{nuccore_to_host[id_]}\n')
