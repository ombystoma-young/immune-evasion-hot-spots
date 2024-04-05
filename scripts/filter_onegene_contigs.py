import argparse

import pandas as pd

from numpy import nan

def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_find_novel_adgs" pipeline. '
                                                 'Remove loci with only one gene')
    parser.add_argument('-g', '--gff', default=None, type=str, nargs='?',
                        help='path to input gff file')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output gff file')
    return parser.parse_args()


# read gff
def _split_attributes(attrs: str):
    pairs = attrs.split(';')
    kvs = [pair.split('=', maxsplit=1) for pair in pairs]
    for i, kv in enumerate(kvs):
        if len(kv) == 1:
            kvs[i - 1][-1] += f'{kv[0]}'
            kvs.remove(kv)
    return {f'ATTRIBUTE_{k}': v for k, v in kvs}


def unites_attributes_col(row: pd.DataFrame, attribute_cols: pd.Index) -> str:
    attributes = []
    for attribute_name in attribute_cols:
        if row[attribute_name] is not nan:
            suffix = '_'.join(attribute_name.split('_')[1:])
            attributes.append(f'{suffix}={row[attribute_name]}')

    return ';'.join(attributes)


def read_gff(gff_path):
    """
    reads gff file into pd dataframe
    :param gff_path: (str) path to gff file
    :return: pd.DataFrame from gff file
    """
    colnames = ['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    entries = []
    with open(gff_path, 'rt') as gff:
        for line in gff:
            entry = line.strip().split('\t')
            entries.append(entry)

    gff_data = pd.DataFrame(entries, columns=colnames)
    gff_data['attribute_dict'] = gff_data['attributes'].apply(_split_attributes)

    norm_attribute = pd.json_normalize(gff_data.attribute_dict)
    gff_data = pd.concat([gff_data, norm_attribute], axis=1)
    # remove temp columns:
    gff_data = gff_data.drop(columns=['attribute_dict', 'attributes'])
    return gff_data


def write_gff(gff: pd.DataFrame, filename: str) -> None:
    """
    writes pandas table into gff file
    :param gff: (pd.DataFrame) with the columns:
    ['chrm', 'start', 'end', 'strand', 'score', 'name', 'attribute_...', 'attribute_...']
    :param filename: (str) name of output file
    :return: None
    """
    gff_df = gff.copy()
    attribute_cols = gff_df.columns[gff_df.columns.str.startswith('ATTRIBUTE_')]
    gff_df['attributes'] = gff_df.apply(lambda row: unites_attributes_col(row, attribute_cols), axis=1)
    gff_df = gff_df[['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']]
    gff_df.to_csv(filename, sep='\t', index=False, header=False)
    print(f'{filename} successfully written')


def filter_contigs(in_gff: pd.DataFrame) -> pd.DataFrame:
    df = in_gff.copy()
    counts = df.groupby('seq_id')['ATTRIBUTE_ID'].agg(['count'])
    has_more_than_one_gene = list(counts.query('count > 1').index.values)
    df = df[df['seq_id'].isin(has_more_than_one_gene)]
    return df


if __name__ == '__main__':
    in_file = parse_args().gff
    out_file = parse_args().output

    gff_df = read_gff(in_file)
    modified_gff_df = filter_contigs(gff_df)
    write_gff(modified_gff_df, out_file)


