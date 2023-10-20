import os

assemblies_dir = 'ncbi_dataset/data'
meta_dir = 'metadata'

db_name = 'phages_genomes_concat'
cluster_prot_dir = 'protein_clusterization'
cluster_prot_by_dataset_dir = os.path.join(cluster_prot_dir, 'datasets')
results_dir = 'results'

find_lost_dir = 'find_lost'

os.makedirs(find_lost_dir, exist_ok=True)

# do not rename this (used in side-scripts):
tdrs_search_dir = 'minimap2_out'
pics_dir = 'pics'
aln_dir = 'tdr_search_aln'
upstream_dir = 'upstream_search'

# create nucleotide db
# 'mmseqs createdb examples/QUERY.fasta queryDB'
# search
# 'mmseqs search queryDB targetDB resultDB tmp'
# convert results
# 'mmseqs convertalis queryDB targetDB resultDB resultDB.m8'

rule all:
    input:
        res_db = os.path.join(find_lost_dir, 'res_db', 'res_db.m8'),
        clp = os.path.join(find_lost_dir,'clp_res','res_db.m8')
        # tsv = os.path.join(find_lost_dir,'res_db','res_db.tsv')

rule extract_cols_bed:
    input:
        bed=os.path.join(upstream_dir,'upstream_fixed.bed'),
    output:
        bed = os.path.join(upstream_dir, 'upstream_fixed_cols.bed')
    shell:
        """
        cat {input} | cut -f 1-6 > {output}
        """


rule extract_upstream_fna:
    input:
        bed = os.path.join(upstream_dir, 'upstream_fixed_cols.bed'),
        fna = os.path.join('blasted', 'phages_genomes_concat.fna')
    output:
        fna = os.path.join(upstream_dir, 'upstreams.fna')
    conda:
        'envs/bedtools.yml'
    shell:
        """
        bedtools getfasta -fi {input.fna} -bed {input.bed} -fo {output.fna}
        """

rule create_mmseqs_nucl_db:
    input:
        fna = os.path.join(upstream_dir, 'upstreams.fna')
    output:
        os.path.join(find_lost_dir, 'query_genomes', 'query_genomes.fnaDB')
    conda:
        'envs/mmseq2.yml'
    shell:
        """
        mmseqs createdb {input} {output}
        """

#
# rule get_coords:
#     input:
#         clu = os.path.join(find_lost_dir, 'res_db', 'res_db.0'),
#         db = os.path.join(find_lost_dir, 'query_genomes', 'query_genomes.fnaDB')
#     output:
#         tsv = os.path.join(find_lost_dir,  'res_db', 'res_db.tsv')
#     params: clu=os.path.join(find_lost_dir,  'res_db', 'res_db')
#     conda:
#         'envs/mmseq2.yml'
#     shell:
#         """
#         mmseqs createtsv {input.db} {input.db} {params.clu} {output.tsv}
#         """

rule create_mmseqs_clp_db:
    input:
        fna = os.path.join(find_lost_dir, 'clp_dom.faa')
    output:
        os.path.join(find_lost_dir, 'query_genomes', 'clp.faaDB')
    conda:
        'envs/mmseq2.yml'
    shell:
        """
        mmseqs createdb {input} {output}
        """

rule search_nucl_against_protein_db:
    input:
        faa_db = os.path.join(find_lost_dir, 'query_genomes', 'query_genomes.fnaDB'),
        clu_db = os.path.join(cluster_prot_dir, 'DB_clu_rep')
    output:
        res_db = os.path.join(find_lost_dir, 'res_db', 'res_db.dbtype')
    params:
        res = os.path.join(find_lost_dir, 'res_db', 'res_db')
    conda:
        'envs/mmseq2.yml'
    shell:
        """
        mmseqs search {input.faa_db} {input.clu_db} {params.res} tmp
        """



rule search_clp_against_protein_db:
    input:
        faa_db = os.path.join(find_lost_dir, 'query_genomes', 'clp.faaDB'),
        clu_db = os.path.join(cluster_prot_dir, 'upstream_proteins', 'upstream_proteins.fnaDB')
    output:
        res_db = os.path.join(find_lost_dir, 'clp_res', 'res_db.dbtype')
    params:
        res = os.path.join(find_lost_dir, 'clp_res', 'res_db')
    conda:
        'envs/mmseq2.yml'
    shell:
        """
        mmseqs search {input.faa_db} {input.clu_db} {params.res} tmp
        """


rule resdb2tsv_clp:
    input:
        faa_db = os.path.join(find_lost_dir, 'query_genomes', 'clp.faaDB'),
        clu_db = os.path.join(cluster_prot_dir, 'upstream_proteins', 'upstream_proteins.fnaDB'),
        res_db = os.path.join(find_lost_dir, 'clp_res', 'res_db.dbtype')
    output:
        res_db = os.path.join(find_lost_dir, 'clp_res', 'res_db.m8')
    params:
        res_db=os.path.join(find_lost_dir, 'clp_res', 'res_db')
    conda:
        'envs/mmseq2.yml'
    shell:
        """
        mmseqs convertalis {input.faa_db} {input.clu_db} {params.res_db} {output}
        """


rule resdb2tsv:
    input:
        fna_db = os.path.join(find_lost_dir, 'query_genomes', 'query_genomes.fnaDB'),
        clu_db = os.path.join(cluster_prot_dir, 'DB_clu_rep'),
        res_db = os.path.join(find_lost_dir, 'res_db', 'res_db.dbtype')
    output:
        res_db = os.path.join(find_lost_dir, 'res_db', 'res_db.m8')
    params:
        res_db=os.path.join(find_lost_dir,'res_db', 'res_db')
    conda:
        'envs/mmseq2.yml'
    shell:
        """
        mmseqs convertalis {input.fna_db} {input.clu_db} {params.res_db} {output}
        """
