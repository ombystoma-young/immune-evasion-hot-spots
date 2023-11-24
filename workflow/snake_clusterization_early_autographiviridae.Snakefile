import os

configfile: 'config_autographiviridae.yaml'
# MMSeqs2 version: 9b9383a3f0c7a99fef37640823174368f74c1487,
# since at conda available version there is a bug
os.makedirs(config['early_prot_db_dir'], exist_ok=True)
os.makedirs(config['early_clu_db_dir'], exist_ok=True)
os.makedirs(config['early_clu_reprs_db_dir'], exist_ok=True)

os.makedirs(config['early_clu_msa_db_dir'], exist_ok=True)
os.makedirs(config['early_hhsuite_db_dir'], exist_ok=True)
os.makedirs(config['early_clu_msa_dir'], exist_ok=True)
os.makedirs(config['early_clans_info_dir'], exist_ok=True)
os.makedirs(config['early_clans_concat_dir'], exist_ok=True)

SMP2 = None

def get_msa_clu_db_names(wildcards):
    # note 1: ck_output is the same as OUTDIR, but this requires
    # the checkpoint to complete before we can figure out what it is!

    # note 2: checkpoints will have attributes for each of the checkpoint
    # rules, accessible by name. Here we use make_some_files
    ck_output = checkpoints.unpack_msa_db.get(**wildcards).output[0]
    SMP, = glob_wildcards(os.path.join(ck_output, "{sample}"))
    return expand(os.path.join(ck_output, "{SAMPLE}"), SAMPLE=SMP)


def get_msa_clu_names(wildcards):
    ck_output = checkpoints.unpack_msa.get(**wildcards).output[0]
    SMP2, = glob_wildcards(os.path.join(ck_output, "{sample}"))
    SMP2.remove('.snakemake_timestamp')
    return expand(os.path.join(ck_output, "{SAMPLE}"), SAMPLE=SMP2)

def get_profiles_names(wildcards):
    ck_output = checkpoints.unpack_msa.get(**wildcards).output[0]
    SMP2, = glob_wildcards(os.path.join(ck_output, "{sample}"))
    SMP2.remove('.snakemake_timestamp')
    return expand(os.path.join(config['early_clans_info_dir'], "{SAMPLE}.tsv"), SAMPLE=SMP2)


rule all:
    input:
        os.path.join(config['early_clu_db_dir'], 'early_proteins_clu.faa'),
        os.path.join(config['early_clu_db_dir'], 'early_proteins_clu.tsv'),
        get_profiles_names,
        os.path.join(config['early_clans_concat_dir'], 'early_clans_info.tsv')



# BLOCK PROTEIN SEQUENCES CLUSTERING
rule create_prot_mmseq_db:
    input:
        os.path.join(config['upstreams_dir'], 'early.faa')
    output:
        os.path.join(config['early_prot_db_dir'], 'early_proteins')
    shell:
        """
        mmseqs createdb {input} {output} 
        """

rule cluster_prot:
    input:
        db = os.path.join(config['early_prot_db_dir'], 'early_proteins')
    output:
        clu = os.path.join(config['early_clu_db_dir'], 'early_proteins_clu')
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
    shell:
        """
        mmseqs createsubdb {input.clu} {input.db} {output}
        """


rule get_upstream_clusters_faa:
    input:
        db = os.path.join(config['early_clu_reprs_db_dir'], 'early_clu_reprs')
    output:
        faa = os.path.join(config['early_clu_db_dir'], 'early_proteins_clu.faa')
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

# BLOCK CLASSIFY PROTEIN CLUSTERS INTO CLANS
## CREATE DB
rule create_msa_filtered:
    input:
        db = os.path.join(config['early_prot_db_dir'], 'early_proteins'),
        db_clu = os.path.join(config['early_clu_db_dir'], 'early_proteins_clu')
    output:
        msa = os.path.join(config['early_clu_msa_db_dir'], 'msa.faDB')
    params:
        minseqs = config['min_seq_msa']
    threads: config['maxthreads']
    shell:
        """
        mmseqs result2msa {input.db} {input.db} {input.db_clu} {output.msa} \
        --msa-format-mode 3 --filter-min-enable {params.minseqs} --threads {threads}
        """

checkpoint unpack_msa_db:
    input:
        msa = os.path.join(config['early_clu_msa_db_dir'],'msa.faDB')
    output:
        directory(config['early_msa_unpacked_db_dir'])
    threads: config['maxthreads']
    shell:
        """
        mmseqs unpackdb {input} {output} --unpack-name-mode 1 --threads {threads}
        """

rule build_ffindex:
    input:
        get_msa_clu_db_names
    output:
        data = temp(os.path.join(config['early_hhsuite_db_dir'], 'msa_msa.ffdata')),
        index = temp(os.path.join(config['early_hhsuite_db_dir'], 'msa_msa.ffindex'))
    conda:
        os.path.join(config['envs'], 'hhsuite.yml')
    params:
            in_dir = os.path.join(config['early_msa_unpacked_db_dir'], '')
    threads: 1
    shell:
        """
        ffindex_build -s {output.data} {output.index} {params.in_dir}
        """

rule create_consensus:
    input:
        data = os.path.join(config['early_hhsuite_db_dir'], 'msa_msa.ffdata'),
        index = os.path.join(config['early_hhsuite_db_dir'], 'msa_msa.ffindex')
    output:
        index = temp(os.path.join(config['early_hhsuite_db_dir'], 'msa_a3m_no_order.ffindex')),
        data = temp(os.path.join(config['early_hhsuite_db_dir'], 'msa_a3m_no_order.ffdata'))
    conda:
        os.path.join(config['envs'], 'hhsuite.yml')
    threads: 1
    shell:
        """
        ffindex_apply {input.data} {input.index} -i {output.index} -d {output.data} -- \
        hhconsensus -M 50 -maxres 65535 -i stdin -oa3m stdout -v 0
        """

rule make_hhm:
    input:
        data = os.path.join(config['early_hhsuite_db_dir'], 'msa_a3m_no_order.ffdata'),
        index = os.path.join(config['early_hhsuite_db_dir'], 'msa_a3m_no_order.ffindex'),
    output:
        data = temp(os.path.join(config['early_hhsuite_db_dir'], 'msa_hhm_no_order.ffdata')),
        index = temp(os.path.join(config['early_hhsuite_db_dir'], 'msa_hhm_no_order.ffindex'))
    conda:
        os.path.join(config['envs'], 'hhsuite.yml')
    threads: 1
    shell:
        """
        ffindex_apply {input.data} {input.index} -i {output.index} -d {output.data} -- hhmake -i stdin -o stdout -v 0
        """

rule make_context_states:
    input:
        data = os.path.join(config['early_hhsuite_db_dir'], 'msa_a3m_no_order.ffdata'),
        index = os.path.join(config['early_hhsuite_db_dir'], 'msa_a3m_no_order.ffindex')
    output:
        data = os.path.join(config['early_hhsuite_db_dir'], 'msa_cs219.ffdata'),
        index = os.path.join(config['early_hhsuite_db_dir'], 'msa_cs219.ffindex')
    params:
        pref_in = os.path.join(config['early_hhsuite_db_dir'], 'msa_a3m_no_order'),
        pref_out = os.path.join(config['early_hhsuite_db_dir'], 'msa_cs219')
    conda:
        os.path.join(config['envs'], 'hhsuite.yml')
    threads: 1
    shell:
        """
        cstranslate -f -x 0.3 -c 4 -I a3m -i {params.pref_in} -o {params.pref_out}
        """

rule extract_order_info:
    input:
        os.path.join(config['early_hhsuite_db_dir'], 'msa_cs219.ffindex')
    output:
        temp(os.path.join(config['early_hhsuite_db_dir'], 'sorting.dat'))
    shell:
        """
        sort -k3 -n -r {input} | cut -f1 > {output}
        """

rule change_order_hhm:
    input:
        dat = os.path.join(config['early_hhsuite_db_dir'], 'sorting.dat'),
        data = os.path.join(config['early_hhsuite_db_dir'], 'msa_{type_}_no_order.ffdata'),
        index = os.path.join(config['early_hhsuite_db_dir'], 'msa_{type_}_no_order.ffindex')
    output:
        data = os.path.join(config['early_hhsuite_db_dir'], 'msa_{type_}.ffdata'),
        index = os.path.join(config['early_hhsuite_db_dir'], 'msa_{type_}.ffindex')
    conda:
        os.path.join(config['envs'], 'hhsuite.yml')
    threads: 1
    shell:
        """
        ffindex_order {input.dat} {input.data} {input.index} {output.data} {output.index}
        """


## CREATE QUERY MMSEQS DB, full
rule create_msa:
    input:
        db = os.path.join(config['early_prot_db_dir'], 'early_proteins'),
        db_clu = os.path.join(config['early_clu_db_dir'], 'early_proteins_clu')
    output:
        msa = os.path.join(config['early_clu_msa_dir'], 'msa.faDB')
    threads: config['maxthreads']
    shell:
        """
        mmseqs result2msa {input.db} {input.db} {input.db_clu} {output.msa} \
        --msa-format-mode 3 --threads {threads}
        """

checkpoint unpack_msa:
    input:
        msa = os.path.join(config['early_clu_msa_dir'], 'msa.faDB')
    output:
        directory(config['early_msa_unpacked_dir'])
    threads: config['maxthreads']
    shell:
        """
        mmseqs unpackdb {input} {output} --unpack-name-mode 1 --threads {threads}
        """

rule run_hhsearch:
    input:
        get_msa_clu_names,
        db = expand(os.path.join(config['early_hhsuite_db_dir'],'msa_{type_}.ffdata'),
               type_=['hhm', 'cs219', 'a3m']),
        a3m = os.path.join(config['early_msa_unpacked_dir'], "{sample}")
    output:
        hhr = os.path.join(config['early_clans_info_dir'], '{sample}.hhr'),
        tab = os.path.join(config['early_clans_info_dir'],'{sample}.tsv')
    params:
        db = os.path.join(config['early_hhsuite_db_dir'],'msa')
    threads: 2
    conda:
        os.path.join(config['envs'],'hhsuite.yml')
    shell:
        """
        hhsearch -cpu {threads} -i {input.a3m} -d {params.db} -o {output.hhr} -blasttab {output.tab}
        """

rule concat_search_results:
    input:
        get_profiles_names
    output:
        os.path.join(config['early_clans_concat_dir'], 'early_clans_info.tsv')
    shell:
        """
        cat {input} > {output}
        """
