import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('in_faa', default=None, nargs='?')
    parser.add_argument('in_list', default=None, nargs='?')
    parser.add_argument('out_faa', default=None, nargs='?')
    parser.add_argument('clu', default=None, nargs='?')
    parser.add_argument('extra_faa', default=None, nargs='?')
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


def read_tsv(tsv_path: str) -> dict:
    nuccore_protein = {}
    with open(tsv_path, 'r') as in_tsv:
        for line in in_tsv:
            entry = line.strip().split('\t')
            chrom = entry[0]
            protein_id = entry[1]
            nuccore_protein[protein_id] = chrom
    return nuccore_protein


def read_fa(protein_ids: dict, faa_path: str):
    sequences = {}
    name = None
    with open(faa_path, 'r') as in_faa:
        for line in in_faa:
            if line.startswith('>'):
                if name in protein_ids.keys():
                    sequences[protein_ids[name]] = seq
                name = line.strip().split(' ')[0][1:]
                seq = ''
            else:
                seq += line.strip()
        if name in protein_ids.keys():
            sequences[protein_ids[name]] = seq
    return sequences


if __name__ == '__main__':
    in_faa = parse_args().in_faa
    in_tsv = parse_args().in_list
    out_faa_ = parse_args().out_faa
    clu = parse_args().clu
    extra_faa = parse_args().extra_faa

    protein_ids = read_tsv(in_tsv)
    seqs = read_fa(protein_ids, in_faa)
    print(len(seqs.keys()))

    with open(out_faa_, 'w') as out_faa:
        for chrom in seqs.keys():
            out_faa.write(f'>{chrom}\n')
            out_faa.write(f'{seqs[chrom]}\n')

    if clu == 'samase':
        with open(extra_faa, 'r') as in_faa:
            with open(out_faa_, 'a') as out_faa:
                for line in in_faa:
                    out_faa.write(line)
