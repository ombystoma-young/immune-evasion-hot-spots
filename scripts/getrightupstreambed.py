import os
import csv
import pandas as pd


def get_genomes(path: str) -> list:
    genomes = []
    with open(path, 'r') as tsv_file:
        for line in tsv_file:
            genome = line.strip().split('\t')[0]
            genomes.append(genome)
    return genomes


def get_top2_largest_intergenics(path: str) -> dict:
    df = pd.read_csv(path, sep='\t', names=['genome', 'start', 'end', 'length'])
    df = df.query('`length` >= 500')
    df = df[df.groupby('genome')['length'].rank(method='first', ascending=False) <= 2]
    df = df.groupby(['genome']).agg(lambda col: list(col))
    intergenics = df.to_dict('index')
    return intergenics


class IntersectionMaster:
    def __init__(self, genomes: list, gff_path: str, tsv_path: str,
                 intergenics: dict, stat_file: str, n_t_a: str, bp=5000):
        """

        :param genomes: (list) of genomes needed to analyse
        :param gff_path: (str) path to representative genomes file gff
        :param tsv_path: (str) path to tsv, containing information about TDRs
        :param stat_file: (str) path to file with length of genomes
        :param intergenics: (dict) dict with top2 intergenic regions
        :param n_t_a: (str) path to nuccore to assembly IDs correspondance, tsv
        :param bp: (int) number of bases as "upstream", if there is no borders, such as TDR or intergenic region
        """
        self.genomes = genomes
        self.gff = gff_path
        self.tsv = tsv_path
        self.intergenics = intergenics
        self.bp = bp
        self.stat_file = stat_file
        self.n_t_a = n_t_a
        self.upstreams = None
        self.rnaps = None
        self.tdrs = None
        self._inter_rnaps = None
        self._last_gene_starts = None
        self._lengths = None
        self.nta = None

    def _read_gff(self):
        self._last_gene_starts = {}
        self.rnaps = {}
        genome_id = None
        start = ''
        with open(self.gff, 'r') as gff_file:
            for line in gff_file:
                cds = line.strip().split('\t')
                if cds[0] in self.genomes:
                    if cds[0] != genome_id:
                        self._last_gene_starts[genome_id] = start
                    genome_id = cds[0]
                    start = cds[3]
                    end = cds[4]
                    strand = cds[6]
                    attribute = cds[-1]
                    if attribute.endswith('RNA polymerase'):
                        attribute_list = attribute.split(';')
                        assembly_id = "_".join(attribute_list[0].split('_')[0:2]).split('=')[1]
                        self.rnaps[genome_id] = {'pol_start': start,
                                                 'pol_end': end,
                                                 'strand': strand,
                                                 'assembly_id': assembly_id,
                                                 'attribute': attribute}
        return self

    def _read_tsv(self):
        self.tdrs = {}
        with open(self.tsv, 'r') as tsv_file:
            for line in tsv_file:
                tdr = line.strip().split('\t')
                genome_id = tdr[0]
                if genome_id in self.genomes:
                    start_1, end_1 = tdr[1], tdr[2]
                    start_2, end_2 = tdr[4], tdr[5]
                    # control for 1+ TDRs in one genome
                    if genome_id in self.tdrs:
                        print(f'Warning! Pay attention to {genome_id} TDRs')
                        prev_start2 = int(self.tdrs[genome_id]['tdr_2_start'])
                        prev_end1 = int(self.tdrs[genome_id]['tdr_1_end'])
                        if int(start_2) - int(end_1) < prev_start2 - prev_end1:
                            continue
                    self.tdrs[genome_id] = {
                        'tdr_1_start': start_1,
                        'tdr_1_end': end_1,
                        'tdr_2_start': start_2,
                        'tdr_2_end': end_2
                    }
        return self.tdrs

    def _get_length(self):
        self._lengths = {}
        with open(self.stat_file, 'r') as stat_f:
            for line in stat_f:
                row = line.strip().split('\t')
                assemby_id = row[1]
                length = row[3]
                self._lengths[assemby_id] = length
        return self

    def __read_nta(self):
        self.nta = {}
        with open(self.n_t_a, 'r') as in_f:
            for line in in_f:
                couple = line.strip().split('\t')
                self.nta[couple[0]] = couple[1]
        return self

    def __manual_curation(self, genome_id):
        # rewrite: take them from downloaded gff
        # outliers = {'NC_047772.1': ['ncbi_dataset/data/GCF_002612265.1/genomic.gff', 'GCF_002612265.1'],
        #             'NC_055824.1': ['ncbi_dataset/data/GCF_009388365.1/genomic.gff', 'GCF_009388365.1'],
        #             'NC_021073.1': ['ncbi_dataset/data/GCF_000906335.1/genomic.gff', 'GCF_000906335.1'],
        #             'MN518894.1': ['ncbi_dataset/data/GCA_013426665.1/genomic.gff', 'GCA_013426665.1'],
        #             'NC_049436.1': ['ncbi_dataset/data/GCF_003308735.1/genomic.gff', 'GCF_003308735.1'],
        #             'NC_049380.1': ['ncbi_dataset/data/GCF_002615685.1/genomic.gff', 'GCF_002615685.1']}
        # if genome_id not in outliers.keys():
        #     raise KeyError(f'No info about path to gff: {genome_id}')
        print(f'Manial curation for {genome_id}')
        assembly_id = self.nta[genome_id]
        folder = list(os.walk(f'ncbi_dataset/data/{assembly_id}'))
        if 'genomic.gff' in folder[0][2]:
            path_to_gff = f'ncbi_dataset/data/{assembly_id}/genomic.gff'
            with open(path_to_gff, 'r') as gff_file:
                for line in gff_file:
                    if not line.startswith('#'):
                        row = line.strip().split('\t')
                        if row[2] == 'CDS':
                            attribute = row[-1].split(';')
                            attribute_dict = {i.split('=')[0]: i.split('=')[1] for i in attribute}
                            if attribute_dict['product'] == 'RNA polymerase':
                                start = row[3]
                                end = row[4]
                                strand = row[6]
                                attribute = row[-1]
                                self.rnaps[genome_id] = {'pol_start': start,
                                                         'pol_end': end,
                                                         'strand': strand,
                                                         'assembly_id': assembly_id,
                                                         'attribute': attribute}
                                return genome_id
            if genome_id not in self.rnaps:
                with open('metadata/list_no_pol', 'a') as out_f:
                    out_f.write('\t'.join([assembly_id, genome_id]))
                    out_f.write('\n')
        else:
            with open('metadata/list_no_anno_no_pol', 'a') as out_f:
                out_f.write('\t'.join([assembly_id, genome_id]))
                out_f.write('\n')

    def _subset_pols_dict(self):
        self._read_gff()
        self._read_tsv()
        self.__read_nta()

        subset_rnaps = set(self.rnaps.keys()) & set(self.tdrs.keys())
        for outlier in set(im.tdrs.keys()) - subset_rnaps:
            if outlier != 'seq_id':
                print(f'Warning! Something wrong with annotation of {outlier}')
                id__ = self.__manual_curation(outlier)
                if id__ is not None:
                    subset_rnaps.add(outlier)
        self._inter_rnaps = {id_: self.rnaps[id_] for id_ in subset_rnaps}

        return subset_rnaps

    def __get_upstream(self, rnap, tdr, g_id):
        upstreams = []
        rnap_start, rnap_end = int(rnap['pol_start']), int(rnap['pol_end'])
        tdrs_1_start, tdrs_1_end = int(tdr['tdr_1_start']), int(tdr['tdr_1_end'])
        tdrs_2_start, tdrs_2_end = int(tdr['tdr_2_start']), int(tdr['tdr_2_end'])
        upstream_strand = rnap['strand']
        assembly_id = rnap['assembly_id']
        length = int(self._lengths[assembly_id])
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

    def __get_upstream_intergenic(self, rnap, nuccore_id):
        upstreams = []
        rnap_start, rnap_end = int(rnap['pol_start']), int(rnap['pol_end'])
        upstream_strand = rnap['strand']
        assembly_id = rnap['assembly_id']
        if nuccore_id in self.intergenics.keys():
            intergenics = self.intergenics[nuccore_id]
        else:
            return None
        length = int(self._lengths[assembly_id])
        if len(intergenics['start']) == 1:
            if upstream_strand == '+':
                distance = rnap_start - intergenics['end'][0]
            else:
                distance = intergenics['start'][0] - rnap_end
            if distance > self.bp + 1000:
                return None
            else:
                intergenic_start = intergenics['start'][0]
                intergenic_end = intergenics['end'][0]
        else:
            if upstream_strand == '+':
                distance_1 = rnap_start - intergenics['end'][0]
                distance_2 = rnap_start - intergenics['end'][1]

                if distance_1 < 0 and distance_2 < 0:
                    distance_1 = length - intergenics['end'][0]
                    distance_2 = length - intergenics['end'][1]
                    if distance_1 < distance_2:
                        if distance_1 > self.bp - rnap_start + 1000:
                            return None
                        intergenic_start = intergenics['start'][0]
                        intergenic_end = intergenics['end'][0]
                    else:
                        if distance_2 > self.bp - rnap_start + 1000:
                            return None
                        intergenic_start = intergenics['start'][1]
                        intergenic_end = intergenics['end'][1]
                elif distance_1 > 0 > distance_2:
                    if distance_1 > self.bp + 1000:
                        return None
                    intergenic_start = intergenics['start'][0]
                    intergenic_end = intergenics['end'][0]
                elif distance_2 > 0 > distance_1:
                    if distance_2 > self.bp + 1000:
                        return None
                    intergenic_start = intergenics['start'][1]
                    intergenic_end = intergenics['end'][1]
                else:
                    if distance_1 < distance_2:
                        if distance_1 > self.bp + 1000:
                            return None
                        else:
                            intergenic_start = intergenics['start'][0]
                            intergenic_end = intergenics['end'][0]
                    else:
                        if distance_2 > self.bp + 1000:
                            return None
                        else:
                            intergenic_start = intergenics['start'][1]
                            intergenic_end = intergenics['end'][1]
            else:
                distance_1 = intergenics['start'][0] - rnap_end
                distance_2 = intergenics['start'][1] - rnap_end
                if distance_1 < 0 and distance_2 < 0:
                    distance_1 = intergenics['start'][0]
                    distance_2 = intergenics['start'][1]
                    if distance_1 < distance_2:
                        if distance_1 > self.bp - rnap_start + 1000:
                            return None
                        intergenic_start = intergenics['start'][0]
                        intergenic_end = intergenics['end'][0]
                    else:
                        if distance_2 > self.bp - rnap_start + 1000:
                            return None
                        intergenic_start = intergenics['start'][1]
                        intergenic_end = intergenics['end'][1]
                elif distance_1 > 0 > distance_2:
                    if distance_1 > self.bp + 1000:
                        return None
                    intergenic_start = intergenics['start'][0]
                    intergenic_end = intergenics['end'][0]
                elif distance_2 > 0 > distance_1:
                    if distance_2 > self.bp + 1000:
                        return None
                    intergenic_start = intergenics['start'][1]
                    intergenic_end = intergenics['end'][1]
                else:
                    if distance_1 < distance_2:
                        if distance_1 > self.bp + 1000:
                            return None
                        else:
                            intergenic_start = intergenics['start'][0]
                            intergenic_end = intergenics['end'][0]
                    else:
                        if distance_2 > self.bp + 1000:
                            return None
                        else:
                            intergenic_start = intergenics['start'][1]
                            intergenic_end = intergenics['end'][1]

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

    def __get_upstream_wo_borders(self, rnap):
        upstreams = []
        rnap_start, rnap_end = int(rnap['pol_start']), int(rnap['pol_end'])
        upstream_strand = rnap['strand']
        assembly_id = rnap['assembly_id']
        length = int(self._lengths[assembly_id])

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
                shift = self.bp - rnap_start
                upstream_2_start = str(length - shift - 1)
                upstream_2_end = str(length)
                upstreams.append(
                    [upstream_2_start, upstream_2_end, '.', '0', upstream_strand, str(rnap_start), str(rnap_end)])
            else:
                raise ValueError("Some shit happened with coordinates")

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
                shift = rnap_end + self.bp - length
                upstream_2_start = str(0)
                upstream_2_end = str(shift)
                upstreams.append(
                    [upstream_2_start, upstream_2_end, '.', '0', upstream_strand, str(rnap_start), str(rnap_end)])
            else:
                raise ValueError("Some shit happened with coordinates")

        upstreams.append([upstream_start, upstream_end, '.', '0', upstream_strand, str(rnap_start), str(rnap_end)])
        return upstreams

    def find_upstreams(self):
        self.upstreams = []
        self._get_length()
        genome_with_tdr_ids = self._subset_pols_dict()
        for genome_id in genome_with_tdr_ids:
            rnap_data = self._inter_rnaps[genome_id]
            tdr_data = self.tdrs[genome_id]
            upstream = self.__get_upstream(rnap_data, tdr_data, genome_id)
            if upstream is not None:
                for u in upstream:
                    u.reverse()
                    u.append(genome_id)
                    u.reverse()
                    self.upstreams.append(u)

        genome_wo_tdr_ids = set(self.rnaps.keys()) - genome_with_tdr_ids
        for genome_id in genome_wo_tdr_ids:
            rnap_data = self.rnaps[genome_id]
            upstream = self.__get_upstream_intergenic(rnap_data, genome_id)
            if upstream is None:
                upstream = self.__get_upstream_wo_borders(rnap_data)
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
    upstream_dir = 'upstream_search'
    tdrs_search_dir = 'minimap2_out'
    stats_dir = 'stats'
    meta_dir = 'metadata'
    intergenics_dir = 'promoters_search'
    integenics_path = os.path.join(intergenics_dir, 'all_intergenic_with_length.tsv')
    gff = os.path.join(upstream_dir, 'representative_genomes.gff')
    tsv = os.path.join(tdrs_search_dir, 'TDRs_modified.tsv')
    stat_file = os.path.join(stats_dir, 'genomes_gc_length.statistics')
    genomes_path = os.path.join(meta_dir, 'genomes_after_curation.tsv')
    nuccore_to_assembly = os.path.join(meta_dir, 'assembly_nuccore.tsv')
    bed = os.path.join(upstream_dir, 'upstream_fixed.bed')

    intergenics = get_top2_largest_intergenics(integenics_path)
    genomes_list = get_genomes(genomes_path)
    # for genome in intergenics.keys():
    #      print(genome, intergenics[genome])
    im = IntersectionMaster(genomes=genomes_list,
                            gff_path=gff, tsv_path=tsv,
                            intergenics=intergenics,
                            stat_file=stat_file,
                            n_t_a=nuccore_to_assembly, bp=10000)
    im.find_upstreams()
    # print(im.rnaps)
    # print(im.upstreams)
    im.write_bed(bed)
