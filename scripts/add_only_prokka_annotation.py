import os
from collections import defaultdict
from collections import Counter


def read_prokka_gff(path: str) -> dict:
    couples = {}
    with open(path, 'r') as d_gff_file:
        for line in d_gff_file:
            entry = line.strip().split('\t')
            if entry[2] == 'CDS':
                our_attribute = {el.split('=')[0]: el.split('=')[1] for el in entry[8].split(';')}
                couples[our_attribute['ID']] = our_attribute['product']
    return couples


def assign_clusters_products(path_clus: str, function: dict, long_file_name: str) -> dict:
    clu_leader_functs_only = defaultdict(list)
    with open(long_file_name, 'w') as long_f:
        with open(path_clus, 'r') as clust_proteins:
            for line in clust_proteins:
                couple = line.strip().split('\t')
                clu_parent = couple[0]
                clu_child = couple[1]
                clu_child_func = function[clu_child]
                long_f.write('\t'.join([clu_parent, clu_child, clu_child_func]))
                long_f.write('\n')
                clu_leader_functs_only[clu_parent].append(clu_child_func)
    return clu_leader_functs_only


def write_wide(clusters: dict, wide_file_name: str) -> None:
    with open(wide_file_name, 'w') as out_f:
        for parent, childs in clusters.items():
            c = Counter(childs)
            counts = str({i: c[i] for i in c})
            out_f.write('\t'.join([parent, str(len(childs)), counts]))
            out_f.write('\n')


def read_faa(faa_path):
    sequences = {}
    name = None
    with open(faa_path, 'r') as in_faa:
        for line in in_faa:
            if line.startswith('>'):
                if name is not None:
                    sequences[name] = seq
                name = line.strip().split(' ')[0][1:]
                seq = ''
            else:
                seq += line.strip()
        sequences[name] = seq
    return sequences


def write_wide_with_sequence(clusters: dict, wide_file_name: str, clu_faa: dict) -> None:
    with open(wide_file_name, 'w') as out_f:
        for parent, childs in clusters.items():
            c = Counter(childs)
            counts = str({i: c[i] for i in c})
            seq = clu_faa[parent]
            out_f.write('\t'.join([parent, str(len(childs)), counts, seq]))
            out_f.write('\n')


if __name__ == '__main__':
    our_gff = os.path.join('../upstream_search', 'representative_genomes.gff')
    clus = os.path.join('../protein_clusterization', 'upstream_proteins_clu.tsv')
    clus_faa = os.path.join('../protein_clusterization', 'upstream_proteins_clu.faa')
    clus_out_long = os.path.join('../results', 'upstream_proteins_clu_long_prokka.tsv')
    clus_out_wide = os.path.join('../results', 'upstream_proteins_clu_wide_prokka.tsv')
    clus_out_wide_seq = os.path.join('../results', 'upstream_proteins_clu_wide_seq_prokka.tsv')

    func_prokka = read_prokka_gff(our_gff)
    faa_dict = read_faa(clus_faa)

    clu_wide = assign_clusters_products(clus, func_prokka, clus_out_long)

    write_wide(clu_wide, clus_out_wide)
    write_wide_with_sequence(clu_wide, clus_out_wide_seq, faa_dict)
