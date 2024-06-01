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
    assembly_to_taxid = {}
    nuccore_to_taxid = {}
    file_path = os.path.join(genomes_dir, 'assembly_data_report.jsonl')
    with open(file_path, 'r') as jsonl:
        for line in jsonl:
            entry = json.loads(line.strip())
            taxid = entry['organism']['taxId']
            assembly_id = entry['accession']
            assembly_to_taxid[assembly_id] = taxid
    for assembly_id in assembly_to_taxid.keys():
        seq_rep = os.path.join(genomes_dir, assembly_id, 'sequence_report.jsonl')
        with open(seq_rep, 'rt') as in_file:
            for line in in_file:
                entry = json.loads(line.strip())
                if assembly_id.startswith('GCA'):
                    nuccore = entry['genbankAccession']
                else:
                    nuccore = entry['refseqAccession']
                nuccore_to_taxid[nuccore] = assembly_to_taxid[assembly_id]
    return nuccore_to_taxid


if __name__ == '__main__':
    gen_dir = parse_args().input_dir
    output = parse_args().output

    nuccore_to_taxid = extract_jsonl(genomes_dir=gen_dir)

    with open(output, 'w') as tsv:
        for id_ in nuccore_to_taxid.keys():
            tsv.write(f'{id_}\t{nuccore_to_taxid[id_]}\n')
