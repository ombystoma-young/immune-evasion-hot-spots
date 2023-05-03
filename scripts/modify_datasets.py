import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('wildcard', default=None, nargs='?')
    return parser.parse_args()


def read_data(in_path: str) -> list:
    genomes = []
    with open(in_path, 'r') as in_file:
        for line in in_file:
            genomes.append(line.strip())
    return genomes

def modify_list(genomes: list, dataset: str) -> list:

    mutants_t7 = ['MZ375237.1',
                    'MZ375240.1',
                    'MZ375257.1',
                    'MZ375260.1',
                    'MZ375267.1',
                    'MZ375268.1',
                    'MZ375271.1',
                    'MZ375279.1',
                    'MZ375282.1',
                    'MZ375288.1',
                    'MZ375289.1',
                    'MZ375290.1',
                    'MZ375292.1',
                    'OL964740.1']
    t7 = 'NC_001604.1'
    if dataset == 'dataset_1':
        genomes_modified = [genome for genome in genomes if genome not in mutants_t7]
    else:
        genomes_modified = genomes.copy()
        genomes_modified.append(t7)
    return genomes_modified


if __name__ == '__main__':
    dataset = parse_args().wildcard
    datasets_dir = 'define_datasets'
    in_file = os.path.join(datasets_dir, f'genomes_{dataset}.txt')
    out_file = os.path.join(datasets_dir, f'{dataset}_genomes_modified.txt')
    genomes_list = read_data(in_file)
    modified_genomes = modify_list(genomes_list, dataset)
    with open(out_file, 'w') as out_file:
        for genome in modified_genomes:
            out_file.write(f'{genome}\n')
