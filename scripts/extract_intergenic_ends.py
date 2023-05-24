import os
from collections import defaultdict


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


def calculate_length(start: int, end: int):
    return end - start


def calculate_lengths(intergenics: dict) -> dict:
    intergenic_with_lengths = defaultdict(list)
    for_deletion = set()
    for chrom in intergenics.keys():
        for intergenic in intergenics[chrom]:
            start = intergenic[0]
            end = intergenic[1]
            length = calculate_length(start, end)
            intergenic_with_lengths[chrom].append((start, end, length))
        if intergenics[chrom][0][0] == 0:
            if len(intergenics[chrom]) == 1:
                for_deletion.add(chrom)
            else:
                intergenic_with_lengths[chrom].pop(0)
                intergenic_with_lengths[chrom].pop(-1)
                start_left = intergenics[chrom][0][0]
                end_left = intergenics[chrom][0][1]
                start_right = intergenics[chrom][-1][0]
                end_right = intergenics[chrom][-1][1]
                length_left = calculate_length(start_left, end_left)
                length_right = calculate_length(start_right, end_right)
                length = length_left + length_right
                intergenic_with_lengths[chrom].append((start_right, end_left, length))
    for defective_genome in for_deletion:
        intergenic_with_lengths.pop(defective_genome)
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
    in_file = os.path.join('promoters_search', 'all_intergenic.bed')
    out_file = os.path.join('promoters_search', 'all_intergenic_with_length.tsv')
    intergens = read_bed(in_file)
    intergens_with_length = calculate_lengths(intergens)
    write_table(intergens_with_length, out_file)