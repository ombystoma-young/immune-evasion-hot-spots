import argparse

from numpy import nan
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_annotate_genomes" pipeline. '
                                                 'Create gff file from phanotate output fasta file (translated), '
                                                 'filter cds based on strandness (optional)')
    parser.add_argument('-f', '--faa', default=None, type=str, nargs='?',
                        help='path to raw fasta')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output gff file')
    parser.add_argument('-t', '--tsv', default=None, type=str, nargs='?',
                        help='path to table with predicted domains')
    parser.add_argument('-c', '--clusters', default=None, type=str, nargs='?',
                        help='path to table with clusters')
    parser.add_argument('-s', '--strandness', action='store_true',
                        help='filter based on strandness')
    return parser.parse_args()


def read_fasta_phanotate(in_path: str) -> pd.DataFrame:
    """
    reads fasta file into dict {seq_id: sequence}
    :param in_path: (str) path to fasta file to read
    :return: (pd.DataFrame) of parsed sequence entries with the following columns:
    'chrm', 'start', 'end', 'strand', 'score', 'name'
    """
    entries = []
    with open(in_path, 'rt') as fasta:
        for i, record in enumerate(SimpleFastaParser(fasta)):
            descr = record[0].split(' ')
            chrm = '.'.join(descr[0].split('.')[:-1])
            rel_end = int(descr[0].split('.')[-1])
            rev_start = int(descr[1].split('=')[1][:-1])
            score = round(float(descr[2].split('=')[1][:-1]), 4)
            strand = '+' if rev_start < rel_end else '-'
            start = min(rev_start, rel_end)
            end = max(rev_start, rel_end)
            name = descr[0]
            entries.append([chrm, start, end, strand, score, name])
    entries_df = pd.DataFrame(entries, columns=['chrm', 'start', 'end', 'strand', 'score', 'attribute_ID'])
    return entries_df


def read_clusters_file(in_path: str) -> dict:
    """
    reads MMSeqs2 clusters tsv, and returns (dict) of correspondence: {child: parent}
    """
    child2parent = {}
    with open(in_path, 'rt') as in_file:
        for line in in_file:
            entry = line.strip().split('\t')
            child2parent[entry[1]] = entry[0]
    return child2parent


def find_strandness_by_counting(gff: pd.DataFrame) -> dict:
    """
    finds strandness of chromosome, counting number of CDSs on + or - strand. Applicable only for directed chromosomes.
    :param gff: (pd.DataFrame)
    :return: (dict): {chromosome: strand}
    """
    counts = gff.groupby(['chrm', 'strand'])['attribute_ID'].count().reset_index()
    idx = counts.groupby('chrm')['attribute_ID'].transform(max) == counts['attribute_ID']
    max_counts = counts[idx][['chrm', 'strand']].set_index('chrm')
    return max_counts.to_dict()


def extract_domains(table_path: str) -> dict:
    colnames = ['query', 'target', 'alnLen', 'seqIdentity', 'eVal',
                'qStart', 'qEnd', 'qLen', 'tStart', 'tEnd', 'tLen']
    df = pd.read_csv(table_path, sep='\t', names=colnames)
    agg_df = df.groupby('target').agg({'query': lambda x: ",".join(list(x))}).reset_index()
    agg_df.rename(columns={"target": "attribute_ID", "query": "PHROGs"}, errors="raise", inplace=True)
    return agg_df.set_index('attribute_ID').to_dict()


def filter_by_strand(gff: pd.DataFrame, strandnesses: dict) -> pd.DataFrame:
    gff_modified = gff.copy()
    gff_modified['pass_strand'] = gff_modified.apply(lambda row: row['strand'] == strandnesses['strand'][row['chrm']], axis=1)
    gff_modified = gff_modified[gff_modified['pass_strand']].drop(columns=['pass_strand'])

    return gff_modified


def find_domain(attribute_id: str, prot2par: dict, par2dom: dict):
    if attribute_id not in prot2par.keys():
        return nan
    elif prot2par[attribute_id] not in par2dom['PHROGs'].keys():
        return nan
    else:
        return par2dom['PHROGs'][prot2par[attribute_id]]


def add_domains_info(gff: pd.DataFrame, prot2par: dict, par2dom: dict) -> pd.DataFrame:
    gff_modified = gff.copy()
    gff_modified['attribute_PHROGs'] = gff_modified.apply(lambda row: find_domain(row['attribute_ID'],
                                                                                  prot2par, par2dom), axis=1)
    return gff_modified


def unites_attributes_col(row: pd.DataFrame, attribute_cols: pd.Index) -> str:
    attributes = []
    for attribute_name in attribute_cols:
        if row[attribute_name] is not nan:
            suffix = '_'.join(attribute_name.split('_')[1:])
            attributes.append(f'{suffix}={row[attribute_name]}')

    return ';'.join(attributes)


def write_gff(gff: pd.DataFrame, filename: str) -> None:
    """
    writes pandas table into gff file
    :param gff: (pd.DataFrame) with the columns:
    ['chrm', 'start', 'end', 'strand', 'score', 'name', 'attribute_...', 'attribute_...']
    :param filename: (str) name of output file
    :return: None
    """
    gff_df = gff.copy()
    attribute_cols = gff_df.columns[gff_df.columns.str.startswith('attribute')]
    gff_df['attributes'] = gff_df.apply(lambda row: unites_attributes_col(row, attribute_cols), axis=1)
    gff_df['type'] = 'CDS'
    gff_df['source'] = 'PHANOTATE'
    gff_df['phase'] = 0
    gff_df = gff_df[['chrm', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']]
    gff_df.to_csv(filename, sep='\t', index=False, header=False)
    print(f'{filename} successfully written')


if __name__ == '__main__':
    faa_path = parse_args().faa
    tsv_path = parse_args().tsv
    clusters = parse_args().clusters
    strandness = parse_args().strandness
    output_path = parse_args().output

    gff = read_fasta_phanotate(faa_path)
    prot2parent = read_clusters_file(in_path=clusters)
    parentid2domains = extract_domains(table_path=tsv_path)

    if strandness:
        strandnesses = find_strandness_by_counting(gff=gff)
        gff = filter_by_strand(gff=gff, strandnesses=strandnesses)

    gff = add_domains_info(gff=gff, prot2par=prot2parent, par2dom=parentid2domains)

    write_gff(gff, output_path)
