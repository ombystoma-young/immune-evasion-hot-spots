import argparse
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_process_clustering_results" pipeline. '
                                                 'Parse CL1 results into tsv file')
    parser.add_argument('-i', '--input', default=None, type=str, nargs='?',
                        help='path to input table file')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output table file')
    return parser.parse_args()


def read_clu1_res(in_path) -> dict:
    clu2clans = defaultdict(list)
    with open(in_path, 'rt') as in_file:
        for i, line in enumerate(in_file):
            clan = line.strip().split('\t')
            name = f'clan_{i+1}'
            for cluster in clan:
                clu2clans[cluster].append(name)
    return clu2clans


def fix_overlaps(clu2clans: dict) -> dict:
    prefix = f'multi_clan_'
    clans = {}
    for cluster in clu2clans.keys():
        if len(clu2clans[cluster]) > 1:
            suffix = '_'.join([clan_name.split('_')[-1] for clan_name in clu2clans[cluster]])
            clans[cluster] = prefix + suffix
        else:
            clans[cluster] = clu2clans[cluster][0]
    return clans


def write_clusters(out_path: str, clu2clans: dict) -> None:
    with open(out_path, 'wt') as out_file:
        for pair in clu2clans.items():
            out_file.write('\t'.join(pair))
            out_file.write('\n')


if __name__ == '__main__':
    in_path = parse_args().input
    out_path = parse_args().output
    clu2clans = read_clu1_res(in_path)
    clu2clans_fixed = fix_overlaps(clu2clans)
    write_clusters(out_path, clu2clans_fixed)
