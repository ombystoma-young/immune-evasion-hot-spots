import argparse
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_filter_early" pipeline. '
                                                 'Find lengths of intergenic regions, '
                                                 'taking into account circular genomes')
    parser.add_argument('-i', '--input', default=None, type=str, nargs='?',
                        help='path to input bed with intergenic regions cooridnates')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output tsv with lengths')
    parser.add_argument('-l', '--lengths', default=None, type=str, nargs='?',
                        help='path to bed file with end eqaul to length of chromosome')
    return parser.parse_args()


def read_chrm_lengths(path_chr_length: str) -> dict:
    chrm_lengths = {}
    with open(path_chr_length, 'rt') as chr_len_file:
        for line in chr_len_file:
            chrm, _, length = line.strip().split('\t')
            chrm_lengths[chrm] = int(length)
    return chrm_lengths


def read_bed(path_bed: str) -> dict:
    intergenics = defaultdict(list)
    with open(path_bed, 'r') as bed_file:
        for line in bed_file:
            entry = line.strip().split('\t')
            chrom = entry[0]
            start = int(entry[1])
            end = int(entry[2])
            intergenics[chrom].append((start, end))
    return intergenics


def is_intergenics_split_by_termini(first_start: int, last_end: int, chrm_len: int):
    if first_start == 0 and last_end == chrm_len:
        return True
    else:
        return False


def calculate_lengths(intergenics: dict, chrm_lengths: dict) -> dict:
    intergenic_with_lengths = defaultdict(list)
    for chrom in intergenics.keys():
        for intergenic in intergenics[chrom]:
            start = intergenic[0]
            end = intergenic[1]
            length = end - start
            intergenic_with_lengths[chrom].append((start, end, length))
        # check if the one big intergenic is split by termini
        first_start = intergenics[chrom][0][0]
        last_end = intergenics[chrom][-1][1]
        chrm_len = chrm_lengths[chrom]
        if is_intergenics_split_by_termini(first_start, last_end, chrm_len):
            # if so, unite two intergenic regions into one big
            intergenic_with_lengths[chrom].pop(0)  # remove the first intergenic
            intergenic_with_lengths[chrom].pop(-1)  # remove the last intergenic
            start_first, end_first = intergenics[chrom][0]
            start_last, end_last = intergenics[chrom][-1]
            length = end_first - start_first + end_last - start_last
            intergenic_with_lengths[chrom].append((start_last, end_first, length))
    return intergenic_with_lengths


def write_table(intergen_with_length: dict, path_out: str) -> None:
    with open(path_out, 'w') as out_file:
        for chrom in intergen_with_length.keys():
            for intergen in intergen_with_length[chrom]:
                start = str(intergen[0])
                end = str(intergen[1])
                length = str(intergen[2])
                string = '\t'.join([chrom, start, end, length])
                out_file.write(string)
                out_file.write('\n')


if __name__ == '__main__':
    in_file = parse_args().input
    out_file = parse_args().output
    lengths_file = parse_args().lengths

    chrm_lengths = read_chrm_lengths(lengths_file)
    intergens = read_bed(in_file)
    intergens_with_length = calculate_lengths(intergens, chrm_lengths)
    write_table(intergens_with_length, out_file)