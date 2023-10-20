
def read_clusters_file(path: str) -> dict:
    clusters2entries = {}
    with open(path, 'rt') as clusters_file:
        for cluster_num, line in enumerate(clusters_file):
            entries = line.strip().split('\t')
            clusters2entries[cluster_num] = entries
    entries2clusters = {}
    for cluster_num, entries in clusters2entries.items():
        for entry in entries:
            entries2clusters[entry] = f'CL_{cluster_num + 1}'
    return entries2clusters


def read_clusters_loci(path: str) -> dict:
    entries2clusters = {}
    with open(path, 'rt') as clusters_file:
        for line in clusters_file:
            entry = line.strip().split(',')
            node = entry[0]
            if entry[2] == 'Clustered':
                cluster = entry[3]
            elif entry[2].startswith('Overlap'):
                cluster = entry[2].split(' ')[1][1:-1]
            elif entry[2] == 'Singleton':
                cluster = 'Singleton'
            else:
                cluster = 'Outlier'
            entries2clusters[node] = cluster
    return entries2clusters


def read_pcs2protid(path: str) -> dict:
    prot_ids2pcs = {}
    with open(path, 'rt') as pc2protid_file:
        for line in pc2protid_file:
            corr = line.strip().split(',')
            if not corr[3] in prot_ids2pcs.keys():
                prot_ids2pcs[corr[3]] = corr[0]
    return prot_ids2pcs


def read_ntwk_file(path: str) -> list:
    edges = []
    with open(path, 'rt') as ntwk_file:
        for line in ntwk_file:
            edge = line.strip().split(' ')
            edges.append(edge)
    return edges


def add_clusters_info(edges: list, entries2clusters: dict) -> list:
    """
    adds to edges info additional columns about clusters for 1st vertex (1st column)
    """
    updated_edges = []
    for edge in edges:
        cluster = entries2clusters[edge[0]]
        updated_edges.append(edge + [cluster])
    return updated_edges


if __name__ == '__main__':
    clu_path = '../try_vconnect/test_output_cl1/modules_mcl_5.0.clusters'
    clu_prot = '../try_vconnect/test_output_cl1/vConTACT_proteins.csv'
    ntwk_path = '../try_vconnect/test_output_cl1/modules.ntwk'
    upd_ntwk_path = '../try_vconnect/test_output_cl1/modules_clusters.ntwk'

    # clu_path = '../try_vconnect/test_output3/genome_by_genome_overview.csv'
    # ntwk_path = '../try_vconnect/test_output3/c1.ntw'
    # upd_ntwk_path = '../try_vconnect/test_output3/c1_upd.ntw'
    clu_prev_long = '../results/upstream_proteins_clu_long.tsv'
    clu_prev_def = '../results/clusters - add tblast and union based on phrOG, no missed.tsv'
    childs2edges = read_pcs2protid(clu_prot)
    clusters = read_clusters_file(clu_path)
    ntwk = read_ntwk_file(ntwk_path)
    ntwk_updated = add_clusters_info(edges=ntwk, entries2clusters=clusters)
    with open(upd_ntwk_path, 'wt') as upd_ntwk_file:
        for upd_edge_info in ntwk_updated:
            upd_ntwk_file.write(' '.join(upd_edge_info))
            upd_ntwk_file.write('\n')