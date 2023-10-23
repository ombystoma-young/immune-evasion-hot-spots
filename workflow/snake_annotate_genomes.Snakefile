import os

configfile: 'config_annotate.yaml'

# CREATE FOLDERS
os.makedirs(config['annotation_dir'], exist_ok=True)
os.makedirs(config['annotation_dir_nucl'], exist_ok=True)
os.makedirs(config['annotation_dir_prot'], exist_ok=True)
os.makedirs(config['prot_db_dir'], exist_ok=True)
os.makedirs(config['domains_tables_dir'], exist_ok=True)

# DEFINE PARAMETERS
if config['genomes_type'] != 'meta':
    meta_mode = ''
else:
    meta_mode = '--meta'

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
        os.path.join(config['domains_tables_dir'], "phrog_whole_proteome.tsv")

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
        os.path.join(config['annotation_dir_prot'], '{sample}.faa')
    threads: 1
    conda: os.path.join(config['envs'], 'biopython.yml')
    params:
        script = os.path.join(config['scripts'], 'filter_and_translate.py'),
        threshold = config['prefilter_ann']
    shell:
        """
        python {params.script} -t {params.threshold} -i {input} -o {output}
        """

rule concatenate_prots:
    input:
        expand(os.path.join(config['annotation_dir_prot'], '{sample}.faa'), sample = samples)
    output:
        os.path.join(config['annotation_dir'], 'concatenated.faa')
    shell:
        """
        cat {input} > {output}
        """


rule create_db_phrog:
    input:
        faa = os.path.join(config['annotation_dir'], 'concatenated.faa')
    output:
        os.path.join(config['prot_db_dir'], 'proteomes')
    conda:
        'envs/mmseq2.yml'
    shell:
        """
        mmseqs createdb {input} {output}
        """

rule search_phrog:
    input:
        db = os.path.join(config['prot_db_dir'], 'proteomes')
    output:
        search = os.path.join(config['prot_db_dir'], 'proteomes_phrog.0')
    conda:
        'envs/mmseq2.yml'
    params:
            out = os.path.join(config['prot_db_dir'], 'upstream_phrog'),
            phrog_db = os.path.join(config['phrog_db'], 'phrogs_mmseqs_db', 'phrogs_profile_db'),
            temp_dir = 'tmp_phrog'
    threads: config['maxthreads']
    shell:
        """
        mmseqs search {params.phrog_db} {input} {params.out} {params.temp_dir} --start-sens 1 --sens-steps 3 -s 7
        """

rule write_results:
    input:
        search = os.path.join(config['prot_db_dir'], 'proteomes_phrog.0'),
        db= os.path.join(config['prot_db_dir'], 'proteomes')
    output:
        os.path.join(config['domains_tables_dir'], "phrog_whole_proteome.tsv")
    params:
            search = os.path.join(config['prot_db_dir'], 'proteomes_phrog'),
            phrog_db = os.path.join(config['phrog_db'],  'phrogs_mmseqs_db', 'phrogs_profile_db')
    conda:
        'envs/mmseq2.yml'
    shell:
        """
        mmseqs createtsv {params.phrog_db} {input.db} {params.search} {output}
        """
