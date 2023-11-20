import argparse
import os

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_process_clustering_results" pipeline. '
                                                 'Read information about clustering by MMSeqs2')
    parser.add_argument('-i', '--input', default=None, type=str, nargs='?',
                        help='path to dir with clusters')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output table file')
    return parser.parse_args()


def read_a3m_file(in_path: str) -> dict:
    clu_n_reprs = {
                   'clu': [],
                   'prot': []
                   }
    clu_num = None
    with open(in_path, 'rt') as in_file:
        for line in in_file:
            if line.startswith('#'):
                clu_num = line.split('|')[0][4:]

            elif line.startswith('>'):
                prot_id = line.strip()[1:]
                clu_n_reprs['prot'].append(prot_id)
                clu_n_reprs['clu'].append(clu_num)
    return clu_n_reprs


def process_msas(in_dir: str) -> pd.DataFrame:
    cluster_reprs = {'clu': [],
                     'prot': []
                     }
    for aln in os.listdir(in_dir):
        if aln != '.snakemake_timestamp':
            path = os.path.join(in_dir, aln)
            one_aln = read_a3m_file(path)
            cluster_reprs['clu'] += one_aln['clu']
            cluster_reprs['prot'] += one_aln['prot']
    df = pd.DataFrame.from_dict(cluster_reprs)
    return df


if __name__ == '__main__':
    inpath = parse_args().input
    outpath = parse_args().output
    df = process_msas(inpath)
    df.to_csv(outpath, sep='\t', header=True, index=False)
