import argparse
import pandas as pd
from numpy import nan


def parse_args():
    parser = argparse.ArgumentParser(description='Part of process clustering results pipeline.'
                                                 'Reformat gff3 entries')
    parser.add_argument('-i', '--input', default=None, type=str, nargs='?',
                        help='path to gff3 fasta')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output gff3 file')
    parser.add_argument('-c', '--clustersinfo', default=None, type=str,
                        help='path to file with information about clusters')
    return parser.parse_args()


def _split_attributes(attrs: str):
    pairs = attrs.split(';')
    kvs = [pair.split('=', maxsplit=1) for pair in pairs]
    for i, kv in enumerate(kvs):
        if len(kv) == 1:
            kvs[i - 1][-1] += f'{kv[0]}'
            kvs.remove(kv)
    return {f'ATTRIBUTE_{k}': v for k, v in kvs}


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
            if not line.startswith('#'):
                entry = line.strip().split('\t')
                entries.append(entry)
    gff_data = pd.DataFrame(entries, columns=colnames)
    gff_data['attribute_dict'] = gff_data['attributes'].apply(_split_attributes)

    norm_attribute = pd.json_normalize(gff_data.attribute_dict)
    gff_data = pd.concat([gff_data, norm_attribute], axis=1)
    # remove temp columns:
    gff_data = gff_data.drop(columns=['attribute_dict', 'attributes'])
    return gff_data


def read_tsv(in_file: str) -> pd.DataFrame:
    df = pd.read_csv(in_file, sep='\t')
    df = df.rename(columns={'prot': 'ATTRIBUTE_ID',
                            'prelim_info': 'ATTRIBUTE_prelim_info',
                            'clu': 'ATTRIBUTE_clu',
                            'clan': 'ATTRIBUTE_clan'
                            })
    df['ATTRIBUTE_clan'] = df['ATTRIBUTE_clan'].str.replace(';', ',')
    return df


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
    ['chrm', 'start', 'end', 'strand', 'score', 'name', 'ATTRIBUTE_...', 'ATTRIBUTE_...']
    :param filename: (str) name of output file
    :return: None
    """
    gff_df = gff.copy()
    attribute_cols = gff_df.columns[gff_df.columns.str.startswith('ATTRIBUTE_')]
    gff_df['attributes'] = gff_df.apply(lambda row: unites_attributes_col(row, attribute_cols), axis=1)
    gff_df = gff_df[['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']]
    gff_df.to_csv(filename, sep='\t', index=False, header=False)
    print(f'{filename} successfully written')


def modify_gff(source_gff_df: pd.DataFrame, new_info_df: pd.DataFrame, key: str) -> pd.DataFrame:
    gff_df = source_gff_df.copy()
    gff_df = gff_df.set_index(key).join(new_info_df.set_index(key)).reset_index()
    return gff_df


if __name__ == '__main__':
    input_path = parse_args().input
    output_path = parse_args().output
    clusters_info_path = parse_args().clustersinfo

    gff_df = read_gff(input_path)
    new_info = read_tsv(clusters_info_path)
    gff_df = modify_gff(gff_df, new_info, 'ATTRIBUTE_ID')
    write_gff(gff_df, output_path)
