import os
from collections import defaultdict


def read_bed(path_bed: str) -> dict:
    upstreams = defaultdict(list)
    with open(path_bed, 'r') as bed_file:
        for line in bed_file:
            entry = line.strip().split('\t')
            if len(entry) == 12:
                chrom = entry[0]
                strand = entry[5]
                start_pol = int(entry[6])
                end_pol = int(entry[7])
                start_tdr_1 = int(entry[8])
                end_tdr_1 = int(entry[9])
                start_tdr_2 = int(entry[10])
                end_tdr_2 = int(entry[11])
                upstream_end = int(entry[2])
                upstreams[chrom].append((strand, start_pol, end_pol,
                                         start_tdr_1, end_tdr_1,
                                         start_tdr_2, end_tdr_2,
                                         upstream_end))
            else:
                break
    return upstreams


def calculate_length_simple(coordinates: list):
    strand = coordinates[0]
    if strand == '+':
        start_pol = coordinates[1]
        start_tdr_left = coordinates[3]
        length = start_pol - start_tdr_left
    else:
        end_pol = coordinates[2]
        end_tdr_right = coordinates[6]
        length = end_tdr_right - end_pol
    return length


def calculate_length_two_parts(coordinates: list):
    strand = coordinates[0][0]
    if strand == '+':
        length_left = coordinates[0][1]
        start_tdr_right = coordinates[0][5]
        length_genome = max(coordinates[0][7], coordinates[1][7])
        length_right = length_genome - start_tdr_right
    else:
        length_left = coordinates[0][4]
        end_pol = coordinates[0][2]
        length_genome = max(coordinates[0][7], coordinates[1][7])
        length_right = length_genome - end_pol

    length = length_left + length_right
    return length


def calculate_lengths(upstreams: dict) -> dict:
    upstream_with_lengths = {}
    for chrom in upstreams.keys():
        if len(upstreams[chrom]) > 1:
            length = calculate_length_two_parts(upstreams[chrom])
        else:
            length = calculate_length_simple(upstreams[chrom][0])
        upstream = list(upstreams[chrom][0][:-1])
        upstream.append(length)
        upstream_with_lengths[chrom] = upstream
    print(upstream_with_lengths)

    return upstream_with_lengths


def write_table(upstream_with_lengths: dict, path_out: str) -> None:
    with open(path_out, 'w') as out_file:
        for chrom in upstream_with_lengths.keys():
            upstream = upstream_with_lengths[chrom]
            strand = upstream[0]
            start_pol = str(upstream[1])
            end_pol = str(upstream[2])
            start_tdr_1 = str(upstream[3])
            end_tdr_1 = str(upstream[4])
            start_tdr_2 = str(upstream[5])
            end_tdr_2 = str(upstream[6])
            length = str(upstream[7])
            string = '\t'.join([chrom, strand, start_pol, end_pol,
                                start_tdr_1, end_tdr_1,
                                start_tdr_2, end_tdr_2, length])
            out_file.write(string)
            out_file.write('\n')


if __name__ == '__main__':
    in_file = os.path.join('upstream_search', 'upstream.bed')
    out_file = os.path.join('upstream_search', 'tdr_pol_dist.tsv')
    tdrs = read_bed(in_file)
    tdrs_pols_dist = calculate_lengths(tdrs)
    write_table(tdrs_pols_dist, out_file)
