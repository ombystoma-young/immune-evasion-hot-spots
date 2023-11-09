import os

configfile: 'config_autographiviridae.yaml'

os.makedirs(config['early_prot_db_dir'], exist_ok=True)
os.makedirs(config['early_clu_db_dir'], exist_ok=True)
os.makedirs(config['early_clu_reprs_db_dir'], exist_ok=True)


rule all:
    input:
        os.path.join(config['early_clu_db_dir'], 'early_proteins_clu.faa'),
        os.path.join(config['early_clu_db_dir'],'early_proteins_clu.tsv')


# # BLOCK PROTEIN SEQUENCES CLUSTERIZATION
#
rule create_prot_mmseq_db:
    input:
        os.path.join(config['upstreams_dir'], 'early.faa')
    output:
        os.path.join(config['early_prot_db_dir'], 'early_proteins')
    conda:
        os.path.join(config['envs'], 'mmseq2.yml')
    shell:
        """
        mmseqs createdb {input} {output} --shuffle
        """

rule cluster_prot:
    input:
        db = os.path.join(config['early_prot_db_dir'], 'early_proteins')
    output:
        clu = os.path.join(config['early_clu_db_dir'], 'early_proteins_clu')
    conda:
        os.path.join(config['envs'], 'mmseq2.yml')
    params:
            maxram = config['maxram'],
            tmp = config['mmseqs_temp'],
            clu_min_cov = config['clu_min_cov'],
            clu_min_seq_id = config['clu_min_seq_id'],
    threads: config['maxthreads']
    shell:
        """
        mmseqs cluster --threads {threads} --max-seqs 300 -k 6 --cluster-mode 0 -s 7 \
        --cov-mode 0 -c {params.clu_min_cov} --min-seq-id {params.clu_min_seq_id} --cluster-reassign \
        --split-memory-limit {params.maxram} {input} {output.clu} {params.tmp}
        """


rule get_prot_clusters_tsv:
    input:
        clu = os.path.join(config['early_clu_db_dir'], 'early_proteins_clu'),
        db = os.path.join(config['early_prot_db_dir'], 'early_proteins')
    output:
        tsv = os.path.join(config['early_clu_db_dir'], 'early_proteins_clu.tsv')
    conda:
        os.path.join(config['envs'], 'mmseq2.yml')
    shell:
        """
        mmseqs createtsv {input.db} {input.db} {input.clu} {output.tsv}
        """

rule create_clu_reprs_subdb:
    input:
        clu = os.path.join(config['early_clu_db_dir'], 'early_proteins_clu'),
        db = os.path.join(config['early_prot_db_dir'], 'early_proteins')
    output:
        db = os.path.join(config['early_clu_reprs_db_dir'], 'early_clu_reprs')
    conda:
        os.path.join(config['envs'], 'mmseq2.yml')
    shell:
        """
        mmseqs createsubdb {input.clu} {input.db} {output}
        """


rule get_upstream_clusters_faa:
    input:
        db = os.path.join(config['early_clu_reprs_db_dir'], 'early_clu_reprs')
    output:
        faa = os.path.join(config['early_clu_db_dir'], 'early_proteins_clu.faa')
    conda:
        os.path.join(config['envs'], 'mmseq2.yml')
    shell:
        """
        mmseqs convert2fasta {input.db} {output.faa}
        """

rule get_frequency_table:
    input:
        os.path.join(config['early_clu_db_dir'], 'early_proteins_clu.tsv')
    output:
        os.path.join(config['early_clu_reprs_db_dir'], 'freq_repres_proteins.txt')
    shell:
        """
        cat {input} | cut -f 1 | sort | uniq -dc | sort -n > {output}
        """