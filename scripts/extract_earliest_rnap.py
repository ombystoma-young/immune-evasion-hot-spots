import argparse


def parse_args():  # checked, okay
    parser = argparse.ArgumentParser(description='Part of "snake_filter_early" pipeline. '
                                                 'Finds upstream loci for selected genomes, '
                                                 'based on information about TDRs and big intergenic regions.')
    parser.add_argument('input', default=None, type=str, nargs='?',
                        help='path to input gff file with RNAPs coordinates')
    parser.add_argument('output', default=None, type=str, nargs='?',
                        help='path to output gff file with RNAPs coordinates, 1 per genome')
    return parser.parse_args()


def filter_gff(in_path: str) -> dict:
    """
    extracts the earliest RNAP coordinates, if there is more than one RNAP in genome
    :param in_path:
    :return: (dict) {chromosome name: rnap coodinates}
    """
    entries_per_chrom = {}
    with open(in_path, 'rt') as gff:
        for line in gff:
            entry = line.strip().split('\t')
            chrm = entry[0]
            if chrm not in entries_per_chrom.keys():
                entries_per_chrom[chrm] = entry
            else:
                strand = entry[6]
                cur_start = int(entry[3])
                cur_end = int(entry[4])
                prev_start = int(entries_per_chrom[chrm][3])
                prev_end = int(entries_per_chrom[chrm][4])
                if strand == '+':
                    if cur_start < prev_start:
                        entries_per_chrom[chrm] = entry
                else:
                    if cur_end > prev_end:
                        entries_per_chrom[chrm] = entry
    return entries_per_chrom


if __name__ == '__main__':
    in_file = parse_args().input
    out_file = parse_args().output
    entries = filter_gff(in_file)
    with open(out_file, 'wt') as out_gff:
        for entry in entries.values():
            out_gff.write('\t'.join(entry))
            out_gff.write('\n')