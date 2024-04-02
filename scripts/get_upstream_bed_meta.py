import csv
import pandas as pd
import argparse

'''
usage:
python get_upstream_bed.py -r {input.rnaps.gff} -t {input.tdrs.tsv} -i {input.intergenics.tsv} \
    -l {input.lengths.tsv} -f {input.filtered.txt} -o {output.bed}
'''


def parse_args():  # checked, okay
    parser = argparse.ArgumentParser(description='Part of "snake_filter_early" pipeline. '
                                                 'Finds upstream loci for selected genomes, '
                                                 'based on information about TDRs and big intergenic regions.')
    parser.add_argument('-r', '--rnapsgff', default=None, type=str, nargs='?',
                        help='path to input gff file with RNAPs coordinates')
    parser.add_argument('-i', '--intergenicstsv', default=None, type=str, nargs='?',
                        help='path to tsv with intergenics coordinates')
    parser.add_argument('-l', '--lengthstsv', default=None, type=str, nargs='?',
                        help='path to bed file with end equal to length of chromosome')
    parser.add_argument('-b', '--bp', default=10000, type=int, nargs='?',
                        help='size of upstream, if there is no intergenics and TDRs')
    # parser.add_argument('-s', '--size', default=500, type=int, nargs='?',
    #                     help='minimal intergenic region size')
    parser.add_argument('-s', '--defshort', default=35000, type=int, nargs='?',
                        help='Length of contig to consider it to be incomplete')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output bed file')
    return parser.parse_args()


def get_genomes(path: str) -> list:  # checked, okay
    genomes = []
    with open(path, 'r') as tsv_file:
        for line in tsv_file:
            genome = line.strip().split('\t')[0]
            genomes.append(genome)
    return genomes


def get_intergenics(path: str, ) -> dict:  # checked, okay
    df = pd.read_csv(path, sep='\t', names=['genome', 'start', 'end', 'length'])
    df = df.groupby(['genome']).agg(lambda col: list(col))
    intergenics = df.to_dict('index')
    return intergenics


class IntersectionMaster:
    def __init__(self, gff_path: str,
                 intergenics: dict,
                 lengths_file: str,
                 shortif: int,
                 bp=5000):
        """
        :param gff_path: (str) path to representative genomes file gff
        :param lengths_file: (str) path to file with length of genomes
        :param intergenics: (dict) dict with top2 intergenic regions
        :param bp: (int) number of bases as "upstream", if there is no borders, such as TDR or intergenic region
        :param shortif: (int) number of bases to consider contig as incomplete
        """
        self.gff = gff_path
        self.intergenics = intergenics
        self.bp = bp
        self.lengths_file = lengths_file
        self._short = short
        self.upstreams = None
        self.rnaps = None
        self.tdrs = None
        self._inter_rnaps = None
        self._lengths = None

    def _read_gff(self):  # checked, okay
        self.rnaps = {}
        with open(self.gff, 'r') as gff_file:
            for line in gff_file:
                cds = line.strip().split('\t')
                genome_id = cds[0]
                start = cds[3]
                end = cds[4]
                strand = cds[6]
                attribute = cds[-1]
                self.rnaps[genome_id] = {'pol_start': start,
                                         'pol_end': end,
                                         'strand': strand,
                                         'attribute': attribute}
        return self

    def _get_lengths(self):  # checked, okay
        self._lengths = {}
        with open(self.lengths_file, 'r') as lengths_f:
            for line in lengths_f:
                row = line.strip().split('\t')
                id_ = row[0]
                length = row[2]
                self._lengths[id_] = length
        return self

    def __get_upstream(self, rnap, tdr, g_id):  # kinda checked, okay
        upstreams = []
        rnap_start, rnap_end = int(rnap['pol_start']), int(rnap['pol_end'])
        tdrs_1_start, tdrs_1_end = int(tdr['tdr_1_start']), int(tdr['tdr_1_end'])
        tdrs_2_start, tdrs_2_end = int(tdr['tdr_2_start']), int(tdr['tdr_2_end'])
        upstream_strand = rnap['strand']
        length = int(self._lengths[g_id])
        # RIP ;)
        if upstream_strand == '+':  # "+"-strand
            if (rnap_start > tdrs_1_end) and (rnap_end < tdrs_2_start):  # RNAP BETWEEN TDRs
                if rnap_start > self.bp:  #
                    upstream_start = str(max(tdrs_1_end, rnap_start - self.bp) - 1)
                    upstream_end = str(rnap_end)
                elif rnap_start == self.bp:
                    upstream_start = str(max(tdrs_1_end, rnap_start - self.bp + 1) - 1)
                    upstream_end = str(rnap_end)
                elif rnap_start < self.bp:
                    upstream_start = str(tdrs_1_end - 1)
                    upstream_end = str(rnap_end)
                else:
                    raise ValueError('+ strand, (rnap_start > tdrs_1_end) and (rnap_end < tdrs_2_start) '
                                     'not enough definition')
            elif (rnap_end < tdrs_1_start) and (rnap_end < tdrs_2_start):  # RNAP BEFORE TDRs
                if rnap_start > self.bp:  # enough space at the left
                    upstream_start = str(rnap_start - self.bp - 1)
                    upstream_end = str(rnap_end)
                elif rnap_start == self.bp:  # enough space at the left, equal
                    upstream_start = str(0)
                    upstream_end = str(rnap_end)
                elif rnap_start < self.bp:  # not enough space at the left, need part from the end
                    upstream_start = str(0)
                    upstream_end = str(rnap_end)
                    shift = self.bp - rnap_start
                    if length - shift >= tdrs_2_end:  # at the end enough space before TDR
                        upstream_start_2 = str(length - shift - 1)
                        upstream_end_2 = str(length)
                    else:  # at the end not enough space before TDR, cut after TDR
                        upstream_start_2 = str(tdrs_2_end - 1)
                        upstream_end_2 = str(length)
                    upstreams.append([upstream_start_2, upstream_end_2, '.', '0',
                                      upstream_strand, str(rnap_start),
                                      str(rnap_end), str(tdrs_1_start),
                                      str(tdrs_1_end), str(tdrs_2_start), str(tdrs_2_end)])
                else:
                    raise ValueError('+ strand, (rnap_end < tdrs_1_start) and (rnap_end < tdrs_2_start) '
                                     'not enough definition')
            elif (rnap_end > tdrs_1_end) and (rnap_end > tdrs_2_end):  # RNAP after TDRs
                upstream_start = str(max(tdrs_2_end, rnap_start - self.bp - 1))
                upstream_end = str(rnap_end)
            else:
                self._inter_rnaps.pop(g_id)
                return None
        elif upstream_strand == '-':  # "-"-strand
            if (rnap_start > tdrs_1_end) and (rnap_end < tdrs_2_start):
                upstream_start = str(rnap_start - 1)
                upstream_end = str(min(rnap_end + self.bp, tdrs_2_start))
            elif (rnap_start < tdrs_1_end) and (rnap_end < tdrs_2_end):
                if rnap_end + self.bp <= length:
                    upstream_start = str(rnap_start - 1)
                    upstream_end = str(rnap_end + self.bp)
                elif rnap_end + self.bp > length:
                    upstream_start = str(rnap_start - 1)
                    upstream_end = str(length)
                    shift = self.bp - (length - rnap_end)
                    if tdrs_1_start >= shift:
                        upstream_start_2 = str(0)
                        upstream_end_2 = str(shift)
                    else:
                        upstream_start_2 = str(0)
                        upstream_end_2 = str(tdrs_1_start)
                    upstreams.append([upstream_start_2, upstream_end_2, '.', '0',
                                      upstream_strand, str(rnap_start),
                                      str(rnap_end), str(tdrs_1_start),
                                      str(tdrs_1_end), str(tdrs_2_start), str(tdrs_2_end)])
                elif (rnap_end < tdrs_1_start) and (rnap_end < tdrs_2_start):
                    upstream_start = str(rnap_start - 1)
                    upstream_end = min(rnap_end + self.bp, tdrs_1_start)
                else:
                    raise ValueError
            else:
                self._inter_rnaps.pop(g_id)
                return None
        else:
            raise ValueError('No info about strand')

        upstreams.append([upstream_start, upstream_end, '.', '0',
                          upstream_strand, str(rnap_start), str(rnap_end),
                          str(tdrs_1_start), str(tdrs_1_end),
                          str(tdrs_2_start), str(tdrs_2_end)])
        return upstreams

    def __get_upstream_intergenic(self, rnap, nuccore_id):  # checked, okay
        """
        finds upstreams in case of presence of intergenic region
        :param rnap:
        :param nuccore_id:
        :return: None, if there is no intergenic or it's far away from RNAP gene,
            list of upstream regions otherwise
        """
        # check if there is an intergenic region
        if nuccore_id in self.intergenics.keys():
            intergenics = self.intergenics[nuccore_id]
        else:
            return None

        upstreams = []
        rnap_start, rnap_end = int(rnap['pol_start']), int(rnap['pol_end'])
        upstream_strand = rnap['strand']
        length = int(self._lengths[nuccore_id])
        if upstream_strand == '+':
            distance = rnap_start - intergenics['end'][0]
        else:
            distance = intergenics['start'][0] - rnap_end
        if distance > self.bp + 1000:
            return None
        else:
            intergenic_start = intergenics['start'][0]
            intergenic_end = intergenics['end'][0]

        if upstream_strand == '+':  # "+"-strand
            if rnap_start > intergenic_end:  # RNAP after intergenic
                upstream_start = str(intergenic_end - 1)
                upstream_end = str(rnap_end)
            elif rnap_start < intergenic_end:  # break between
                upstream_start = str(0)
                upstream_end = str(rnap_end)
                upstream_start_2 = str(intergenic_end - 1)
                upstream_end_2 = str(length)
                upstreams.append([upstream_start_2, upstream_end_2, '.', '0',
                                  upstream_strand, str(rnap_start), str(rnap_end),
                                  str(intergenic_start), str(intergenic_end)])
            else:
                raise ValueError(f'{nuccore_id}: problems with intergenics')
        else:  # "-"-strand
            if rnap_end < intergenic_start:  # RNAP after intergenic
                upstream_start = str(rnap_start - 1)
                upstream_end = str(intergenic_start)
            elif rnap_start > intergenic_end:  # break between
                upstream_start = str(rnap_start - 1)
                upstream_end = str(min(length, rnap_end + self.bp))
                if length <= rnap_end+self.bp:
                    upstream_start_2 = str(0)
                    upstream_end_2 = str(intergenic_end)
                    upstreams.append([upstream_start_2, upstream_end_2, '.', '0',
                                      upstream_strand, str(rnap_start), str(rnap_end),
                                      str(intergenic_start), str(intergenic_end)])
            else:
                raise ValueError(f'{nuccore_id}: problems with intergenics')
        upstreams.append([upstream_start, upstream_end, '.', '0',
                          upstream_strand, str(rnap_start), str(rnap_end),
                          str(intergenic_start), str(intergenic_end)])
        return upstreams

    def __get_upstream_wo_borders(self, rnap, nuccore_id):  # checked, okay
        """
        finds upstream of RNAP gene, if there is no borders (TDRs, intergenic region) specified
        :param rnap (dict): coordinates of RNAP ('pol_start', 'pol_end', 'strand', 'attribute')
        :param nuccore_id (str): genome id to link RNAP gene and length of the chromosome
        :return: (list): of upstream regions (two if upstream is split by termini, otherwise one)
        """
        upstreams = []
        rnap_start, rnap_end = int(rnap['pol_start']), int(rnap['pol_end'])
        upstream_strand = rnap['strand']
        length = int(self._lengths[nuccore_id])
        if upstream_strand == '+':  # + strand
            if rnap_start > self.bp:
                upstream_start = str(rnap_start - self.bp - 1)
                upstream_end = str(rnap_end)
            elif rnap_start == self.bp:
                upstream_start = str(rnap_start - self.bp)
                upstream_end = str(rnap_end)
            elif rnap_start < self.bp:
                upstream_start = str(0)
                upstream_end = str(rnap_end)
                if length >= self._short:
                    shift = self.bp - rnap_start
                    upstream_2_start = str(length - shift - 1)
                    upstream_2_end = str(length)
                    upstreams.append(
                        [upstream_2_start, upstream_2_end, '.', '0', upstream_strand, str(rnap_start), str(rnap_end)])
            else:
                raise ValueError("Something happened with coordinates")

        else:  # - strand
            if length > rnap_end + self.bp:
                upstream_start = str(rnap_start - 1)
                upstream_end = str(rnap_end + self.bp)
            elif length == rnap_end + self.bp:
                upstream_start = str(rnap_start - 1)
                upstream_end = str(length)
            elif length < rnap_end + self.bp:
                upstream_start = str(rnap_start - 1)
                upstream_end = str(length)
                if length > self._short:
                    shift = rnap_end + self.bp - length
                    upstream_2_start = str(0)
                    upstream_2_end = str(shift)
                    upstreams.append(
                        [upstream_2_start, upstream_2_end, '.', '0', upstream_strand, str(rnap_start), str(rnap_end)])
            else:
                raise ValueError("Something happened with coordinates")
        upstreams.append([upstream_start, upstream_end, '.', '0', upstream_strand, str(rnap_start), str(rnap_end)])
        return upstreams

    def find_upstreams(self):
        self.upstreams = []
        # determine lengths of chromosomes
        self._get_lengths()
        self._read_gff()

        for genome_id in self.rnaps.keys():
            rnap_data = self.rnaps[genome_id]
            upstream = self.__get_upstream_intergenic(rnap_data, genome_id)
            if upstream is None:
                upstream = self.__get_upstream_wo_borders(rnap_data, genome_id)
            for u in upstream:
                u.reverse()
                u.append(genome_id)
                u.reverse()
                self.upstreams.append(u)

        return self.upstreams

    def write_bed(self, file_name):
        with open(file_name, 'w') as out_f:
            bed_writer = csv.writer(out_f, delimiter='\t')
            for row in self.upstreams:
                bed_writer.writerow(row)
        return self


if __name__ == '__main__':

    gff = parse_args().rnapsgff
    integenics_path = parse_args().intergenicstsv
    lengths_file = parse_args().lengthstsv
    dist = parse_args().bp
    bed = parse_args().output
    short = parse_args().derfshort

    intergenics = get_intergenics(integenics_path)
    im = IntersectionMaster(gff_path=gff,
                            intergenics=intergenics,
                            lengths_file=lengths_file, bp=dist,
                            shortif=short)
    im.find_upstreams()
    im.write_bed(bed)
