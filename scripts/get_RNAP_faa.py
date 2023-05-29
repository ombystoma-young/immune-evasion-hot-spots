import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('wildcard', default=None, nargs='?')
    return parser.parse_args()


def define_genomes_in_dataset(txt_path: str, all=False) -> set:
    genomes = set()
    with open(txt_path, 'r') as in_f:
        for line in in_f:
            if all:
                entry = line.strip().split('\t')[0]
            else:
                entry = line.strip()
            genomes.add(entry)
    return genomes


def read_bed(tsv_path: str, genomes_in_dataset: set) -> dict:
    pol_coordinates = {}
    with open(tsv_path, 'r') as in_tsv:
        for line in in_tsv:
            entry = line.strip().split('\t')
            chrom = entry[0]
            if chrom in genomes_in_dataset:
                pol_start = entry[2]
                pol_end = entry[3]
                pol_coordinates[chrom] = (pol_start, pol_end)
    return pol_coordinates


def extract_from_gff(gff_path: str, pol_coords: dict) -> dict:
    pol_ids = {}
    got = set()
    with open(gff_path, 'r') as in_gff:
        for line in in_gff:
            entry = line.strip().split()
            if entry[0] in pol_coords.keys():
                if entry[0] not in got:
                    chrom = entry[0]
                    start = entry[3]
                    end = entry[4]
                    if start == pol_coords[chrom][0] and end == pol_coords[chrom][1]:
                        got.add(chrom)
                        pol_id = entry[8].split(';')[0].split('=')[1]
                        pol_ids[pol_id] = chrom
    return pol_ids


def read_fa(pol_ids: dict, faa_path: str):
    sequences = {}
    name = None
    with open(faa_path, 'r') as in_faa:
        for line in in_faa:
            if line.startswith('>'):
                if name in pol_ids.keys():
                    sequences[pol_ids[name]] = seq
                name = line.strip().split(' ')[0][1:]
                seq = ''
            else:
                seq += line.strip()
        if name in pol_ids.keys():
            sequences[pol_ids[name]] = seq
    return sequences


if __name__ == '__main__':
    dataset = parse_args().wildcard
    in_tsv = os.path.join('define_datasets', 'joined.tsv')
    in_gff = os.path.join('upstream_search', 'representative_genomes.gff')
    in_faa = os.path.join('upstream_search', 'all_genomes.faa')
    #dataset = 'dataset_3'
    if dataset != 'all':
        in_datasets = os.path.join('define_datasets', f'{dataset}_genomes_modified.txt')
        genomes = define_genomes_in_dataset(in_datasets)
    else:
        in_datasets = os.path.join('metadata', 'genomes_after_curation.tsv')
        genomes = define_genomes_in_dataset(in_datasets, all=True)
    print(genomes)
    pol_coords = read_bed(in_tsv, genomes)
    print(len(pol_coords))
    pol_ids = extract_from_gff(in_gff, pol_coords)
    seqs = read_fa(pol_ids, in_faa)
    print(len(seqs.keys()))
    out_faa = os.path.join('define_datasets', 'alignments', f'polymerases_{dataset}.faa')
    with open(out_faa, 'w') as out_faa:
        for chrom in seqs.keys():
            out_faa.write(f'>{chrom}\n')
            out_faa.write(f'{seqs[chrom]}\n')