import argparse

import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_taxids" pipeline. '
                                                 'Concatenate taxonkit parsed results with nuccore IDs')
    parser.add_argument('-i', '--input_lineage', default='../metadata/lineage/autographiviridae_phage_parsed.lineage',
                        type=str, nargs='?',
                        help='path to taxonkit lineage results',
                        )
    parser.add_argument('-t', '--table', default='../metadata/lineage/autographiviridae_phage.tax.ids',
                        type=str, nargs='?',
                        help='')
    parser.add_argument('-o', '--output', default='../metadata/lineage/autographiviridae_phage_parsed_joined_lineage.tsv',
                        type=str, nargs='?',
                        help='path to output tsv file')
    parser.add_argument('-s', '--hosts', default=None, type=str, nargs='?',
                        help='path to table for hosts')
    return parser.parse_args()


if __name__ == '__main__':
    input_file = parse_args().input_lineage
    tab_file = parse_args().table
    out_file = parse_args().output
    hosts_file = parse_args().hosts

    lineage_df = pd.read_csv(input_file, sep='\t', dtype='str')

    if 'phage' in input_file:
        col_names = ['nuccore_id', 'taxid']
        tab_df = pd.read_csv(tab_file, names=col_names, sep='\t', dtype='str')
        tab_df =tab_df.merge(lineage_df, how='left', on='taxid')
        tab_df.to_csv(out_file, sep='\t', index=False)
    else:
        col_names_tax = ['host_name', 'taxid']
        col_names_nuc = ['nuccore_id', 'host_name']

        tax_tab_df = pd.read_csv(tab_file, names=col_names_tax, sep='\t', dtype='object')
        nuc_tab_df = pd.read_csv(hosts_file, names=col_names_nuc, sep='\t', dtype='object')
        lineage_df = lineage_df.query('superkingdom != "Eukaryota"')
        tax_tab_df = tax_tab_df.merge(lineage_df, how='left', on='taxid')
        nuc_tab_df.merge(tax_tab_df, how='left', on='host_name').drop_duplicates().to_csv(out_file, sep='\t', index=False)


