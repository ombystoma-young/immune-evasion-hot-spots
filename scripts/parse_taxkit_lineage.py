import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_taxids" pipeline. '
                                                 'Parses taxonkit output (-R specified) into table')
    parser.add_argument('-i', '--input_lineage', default=None, type=str, nargs='?',
                        help='path to taxonkit lineage results')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output tsv file')
    parser.add_argument('-m', '--mode', default=None, type=str, nargs='?',
                        help='organism mode: bac, vir')
    return parser.parse_args()


def parse_tk_lineage(in_f: str, out_f: str) -> None:
    entries = {}
    ranks_total = []
    with open(in_f, 'r') as lineage:
        for line in lineage:
            entry = line.strip().split('\t')
            if len(entry) == 1:
                continue
            tax_id = entry[0]
            ranks = entry[2].split(';')
            rank_names = entry[1].split(';')
            entries[tax_id] = {}
            for i, rank in enumerate(ranks):
                if rank not in ranks_total:
                    ranks_total.append(rank)
                entries[tax_id][rank] = rank_names[i]
    with open(out_f, 'w') as tsv:
        tsv.write('taxid\t')
        tsv.write('\t'.join(ranks_total))
        tsv.write('\n')
        for tax in entries.keys():
            tsv.write(tax)
            tsv.write('\t')
            for rank in ranks_total:
                if rank in entries[tax].keys():
                    tsv.write(entries[tax][rank])
                else:
                    tsv.write('NA')
                tsv.write('\t')
            tsv.write('\n')


if __name__ == '__main__':
    in_f = parse_args().input_lineage
    out_f = parse_args().output
    mode = parse_args().mode
    parse_tk_lineage(in_f, out_f)
