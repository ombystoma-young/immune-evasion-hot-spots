import os

configfile: 'config_autographiviridae_refseq.yaml'

# CREATE FOLDERS
os.makedirs(config['annotation_dir'], exist_ok=True)
os.makedirs(config['annotation_dir_nucl'], exist_ok=True)
os.makedirs(config['annotation_dir_prot'], exist_ok=True)
os.makedirs(config['prot_db_dir'], exist_ok=True)
os.makedirs(config['domains_tables_dir'], exist_ok=True)
os.makedirs(config['mmseqs_db'], exist_ok=True)
os.makedirs(config['mmseqs_temp'], exist_ok=True)


# DEFINE PARAMETERS
if config['genomes_type'] != 'meta':
    meta_mode = ''
else:
    meta_mode = '--meta'

if config['filter_strand']:
    strandness = '-s'
else:
    strandness = ''

# READ SAMPLES
samples = []
for sample in os.listdir(config['genomes_source']):
    if config['genomes_type'] == 'refseq':
        if sample.startswith('GC'):
            samples.append(sample)
    else:
        samples.append(sample)

rule all:
    input:
        os.path.join(config['domains_tables_dir'], "phrog_filtered_proteome.tsv"),
        os.path.join(config['prot_clu_db_dir'], 'proteomes_clu.tsv'),
        os.path.join(config['annotation_dir'], 'concatenated.gff')

rule predict_orfs:
    input:
        jsonl = os.path.join(config['genomes_source'], '{sample}', 'sequence_report.jsonl')
    output:
        os.path.join(config['annotation_dir_nucl'], '{sample}.fna')
    params:
        in_fna = lambda wc: os.path.join(config['genomes_source'], f'{wc.sample}')
    threads: 1
    conda:
        os.path.join(config['envs'], 'phanotate.yml')
    shell:
        """
        phanotate.py {params.in_fna}/*.fna -f fasta -o {output}
        """

rule translate_orfs:
    input:
        os.path.join(config['annotation_dir_nucl'], '{sample}.fna')
    output:
        nof = os.path.join(config['annotation_dir_prot'], '{sample}.faa'),
        f = os.path.join(config['annotation_dir_prot'],'filtered_{sample}.faa')
    threads: 1
    conda: os.path.join(config['envs'], 'biopython.yml')
    params:
        script = os.path.join(config['scripts'], 'filter_and_translate.py'),
        threshold = config['prefilter_ann']
    shell:
        """
        python {params.script} -t 0 -i {input} -o {output.nof}
        python {params.script} -t {params.threshold} -i {input} -o {output.f}
        """

rule concatenate_prots:
    input:
        filtered = expand(os.path.join(config['annotation_dir_prot'], 'filtered_{sample}.faa'), sample = samples),
        not_filtered = expand(os.path.join(config['annotation_dir_prot'], '{sample}.faa'), sample = samples)
    output:
        filtered = os.path.join(config['annotation_dir'], 'filtered_concatenated.faa'),
        not_filtered = os.path.join(config['annotation_dir'], 'concatenated.faa'),
    shell:
        """
        cat {input.filtered} > {output.filtered}
        cat {input.not_filtered} > {output.not_filtered}
        """


rule create_db_whole_proteome:
    input:
        faa = os.path.join(config['annotation_dir'], 'filtered_concatenated.faa')
    output:
        os.path.join(config['prot_db_dir'], 'proteomes_filtered')
    conda:
         os.path.join(config['envs'], 'mmseq2.yml')
    shell:
        """
        mmseqs createdb {input} {output}
        """

rule cluster_db:
    input:
        os.path.join(config['prot_db_dir'], 'proteomes_filtered')
    output:
        os.path.join(config['prot_clu_db_dir'], 'proteomes_clu.0')
    params:
        out =  os.path.join(config['prot_clu_db_dir'], 'proteomes_clu'),
        tmp = os.path.join(config['mmseqs_temp'], 'clu_temp'),
        maxram = config['maxram'],
        maxmem = config['maxmem'],
    threads: config['maxthreads']
    conda: os.path.join(config['envs'], 'mmseq2.yml')
    shell:
        """
        mmseqs cluster  \
        --threads {threads} --max-seqs 300 -k 6 --cluster-mode 1 \
        --cov-mode 0 -c 0.7 \
        --split-memory-limit {params.maxram} {input} {params.out} {params.tmp}
        """

rule extract_reprs_db:
    input:
        db = os.path.join(config['prot_db_dir'], 'proteomes_filtered'),
        db_clu = os.path.join(config['prot_clu_db_dir'], 'proteomes_clu.0')
    output:
        tsv = os.path.join(config['prot_clu_db_dir'], 'proteomes_clu.tsv')
    params:
        db_clu = os.path.join(config['prot_clu_db_dir'],'proteomes_clu')
    conda: os.path.join(config['envs'], 'mmseq2.yml')
    shell:
        """
        mmseqs createtsv {input.db} {input.db} {params.db_clu} {output.tsv}
        """

rule create_reprs2search:
    input:
        db = os.path.join(config['prot_db_dir'], 'proteomes_filtered'),
        db_clu = os.path.join(config['prot_clu_db_dir'], 'proteomes_clu.0')
    output:
        db_rep = os.path.join(config['prot_clu_reprs_db_dir'], 'proteome_reprs')
    params:
        db_clu = os.path.join(config['prot_clu_db_dir'], 'proteomes_clu')
    conda: os.path.join(config['envs'],'mmseq2.yml')
    shell:
        """
        mmseqs createsubdb {params.db_clu} {input.db} {output.db_rep}
        """


rule search_phrog:
    input:
        db = os.path.join(config['prot_clu_reprs_db_dir'], 'proteome_reprs')
    output:
        search = os.path.join(config['prot_clu_reprs_db_dir'], 'proteomes_phrog')
    conda:
        os.path.join(config['envs'], 'mmseq2.yml')
    params:
            phrog_db = os.path.join(config['phrog_db'], 'phrogs_mmseqs_db', 'phrogs_profile_db'),
            temp_dir = os.path.join(config['mmseqs_temp'], 'tmp_phrog'),
            maxmem = config['maxmem'],
            maxram = config['maxram']
    threads: config['maxthreads']
    shell:
        """
        mmseqs search \
        --split-memory-limit {params.maxram} --disk-space-limit {params.maxmem} --threads {threads} \
        --start-sens 1 --sens-steps 3 -s 7 \
        {params.phrog_db} {input} {output.search} {params.temp_dir}
        """

rule write_results:
    input:
        search = os.path.join(config['prot_clu_reprs_db_dir'], 'proteomes_phrog'),
        db = os.path.join(config['prot_clu_reprs_db_dir'], 'proteome_reprs')
    output:
        os.path.join(config['domains_tables_dir'], "phrog_filtered_proteome.tsv")
    params:
            search = os.path.join(config['prot_clu_reprs_db_dir'], 'proteomes_phrog'),
            phrog_db = os.path.join(config['phrog_db'],  'phrogs_mmseqs_db', 'phrogs_profile_db')
    conda:
        os.path.join(config['envs'], 'mmseq2.yml')
    shell:
        """
        mmseqs createtsv {params.phrog_db} {input.db} {params.search} {output}
        """

rule faa_to_gff:
    input:
        faa = os.path.join(config['annotation_dir_prot'], '{sample}.faa'),
        domains_tsv = os.path.join(config['domains_tables_dir'], "phrog_filtered_proteome.tsv"),
        clu_parents_tsv = os.path.join(config['prot_clu_db_dir'], 'proteomes_clu.tsv')

    output:
        os.path.join(config['annotation_dir_table'], '{sample}.gff')
    conda:
        os.path.join(config['envs'], 'biopython.yml')
    params:
        script = os.path.join(config['scripts'], 'create_proteome_gff.py'),
        s = strandness
    shell:
        """
        python {params.script} -f {input.faa} -t {input.domains_tsv} -c {input.clu_parents_tsv} -o {output} {params.s}
        """

rule unite_gff:
    input:
        expand(os.path.join(config['annotation_dir_table'], '{sample}.gff'), sample=samples)
    output:
        os.path.join(config['annotation_dir'], 'concatenated.gff')
    shell:
        """
        cat {input} > {output}
        """
