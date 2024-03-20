import os

configfile: 'config_autographiviridae_meta.yaml'

# CREATE FOLDERS
os.makedirs(config['data'], exist_ok=True)
os.makedirs(config['contigs_dir'], exist_ok=True)
os.makedirs(config['checkv_output'], exist_ok=True)
os.makedirs(config['annotation_dir'], exist_ok=True)
os.makedirs(config['annotation_dir_nucl'], exist_ok=True)
os.makedirs(config['annotation_dir_prot'], exist_ok=True)
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

databases = []
paths = {}
actions = {}
with open(os.path.join(config['meta'], config['samples_description']), 'rt') as samples_f:
    for line in samples_f:
        if not line.startswith('#'):
            entry = line.strip().split('\t')
            db_name = entry[0]
            paths[db_name] =  entry[1]
            actions[db_name] = entry[2]
            databases.append(db_name)


def get_faa_files_names(wildcards):
    # note 1: ck_output is the same as OUTDIR, but this requires
    # the checkpoint to complete before we can figure out what it is!

    # note 2: checkpoints will have attributes for each of the checkpoint
    # rules, accessible by name. Here we use make_some_files
    ck_output = checkpoints.move_all_contigs_in_one_folder.get(**wildcards).output[0]
    SMP, = glob_wildcards(os.path.join(ck_output, "{sample}.fna"))
    return expand(os.path.join(config['annotation_dir_prot'], "{SAMPLE}.faa"), SAMPLE=SMP)


def get_faa_filtered_files_names(wildcards):
    # note 1: ck_output is the same as OUTDIR, but this requires
    # the checkpoint to complete before we can figure out what it is!

    # note 2: checkpoints will have attributes for each of the checkpoint
    # rules, accessible by name. Here we use make_some_files
    ck_output = checkpoints.move_all_contigs_in_one_folder.get(**wildcards).output[0]
    SMP, = glob_wildcards(os.path.join(ck_output, "{sample}.fna"))
    return expand(os.path.join(config['annotation_dir_prot'], "filtered_{SAMPLE}.faa"), SAMPLE=SMP)


rule all:
    input:
        expand(os.path.join(config['contigs_dir'],'{db}'), db=databases),
        os.path.join(config['annotation_dir'], 'concatenated.faa'),
        os.path.join(config['annotation_dir'], 'concatenated.gff'),
        config['genomes_source']

# rule setup_db:
#     output:
#         os.path.join(config['checkv_db_path'], 'checkv-db-v1.5', 'README.txt')
#     params:
#         config['checkv_db_path']
#     conda:
#         os.path.join(config['envs'], 'checkv.yml')
#     shell:
#         """
#         checkv download_database {params}
#         """

# rule check_completeness:
#     input:
#         sequences = lambda wc: os.path.join(config['db_fastas'], f'{wc.db}', f'{paths[wc.db]}'),
#         db = os.path.join(config['checkv_db_path'], 'checkv-db-v1.5', 'README.txt')
#     output:
#         os.path.join(config['checkv_output'], '{db}', 'completeness.tsv')
#     conda:
#         os.path.join(config['envs'], 'checkv.yml')
#     params:
#         out_path = os.path.join(config['checkv_output'], '{db}'),
#         db_path = os.path.join(config['checkv_db_path'], 'checkv-db-v1.5')
#     threads: config['maxthreads']
#     shell:
#         """
#         checkv completeness {input.sequences} {params.out_path} -t {threads} -d {params.db_path}
#         """

rule parse_len:
    input:
        os.path.join(config['checkv_output'], '{db}', 'completeness.tsv')
    output:
        os.path.join(config['checkv_output'], '{db}', 'completeness_filtered.tsv')
    params: minlen = config['min_contig_len']
    shell:
        """
        cat {input} | awk -v FS="\t" -v OFS="\t" -v l={params.minlen} '$2>l {{print $0}}' > {output}
        """


rule parse_completeness_results:
    input:
        os.path.join(config['checkv_output'], '{db}', 'completeness_filtered.tsv')
    output:
        os.path.join(config['checkv_output'], '{db}', 'complete_contigs.txt')
    params:
        script = os.path.join(config['scripts'], 'parse_checkv_completeness.py'),
        thres = config['checkv_completeness_thres']
    shell:
        """
        python {params.script} -i {input} -o {output} -t {params.thres}
        """

rule split_fasta:
    input:
        sequences = lambda wc: os.path.join(config['db_fastas'], f'{wc.db}', f'{paths[wc.db]}')
    output:
        directory(os.path.join(config['contigs_dir'], '{db}'))
    params:
        script = os.path.join(config['scripts'], 'split_fastas_meta.py'),
        batches = int(config['maxthreads']) * 5,
        minlen = config['min_contig_len']
    conda: os.path.join(config['envs'], 'biopython.yml')
    shell:
        """
        python {params.script} -i {input.sequences} -o {output} -b {params.batches} -l {params.minlen}
        """

checkpoint move_all_contigs_in_one_folder:
    input:
        expand(os.path.join(config['contigs_dir'], '{db}'), db=databases)
    output:
        directory(config['all_contigs_dir'])
    shell:
        """
        mkdir {output};
        array=({input})
        for dir in ${{array[*]}}
          do cp ${{dir}}/* {output}/
        done
        """

rule predict_orfs:
    input:
        os.path.join(config['all_contigs_dir'], '{sample}.fna')
    output:
        os.path.join(config['annotation_dir_nucl'], '{sample}.fna')
    threads: 1
    conda:
        os.path.join(config['envs'], 'phanotate.yml')
    shell:
        """
        phanotate.py {input} -f fasta -o {output}
        """

rule translate_orfs:
    input:
        rules.predict_orfs.output[0] # os.path.join(config['annotation_dir_nucl'], '{sample}.fna')
    output:
        nof = os.path.join(config['annotation_dir_prot'], '{sample}.faa'),
        f = os.path.join(config['annotation_dir_prot'], 'filtered_{sample}.faa')
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
        filtered = get_faa_filtered_files_names,
        not_filtered = get_faa_files_names
    output:
        filtered = os.path.join(config['annotation_dir'], 'filtered_concatenated.faa'),
        not_filtered = os.path.join(config['annotation_dir'], 'concatenated.faa'),
    shell:
        """
        cat {input.filtered} >> {output.filtered}
        cat {input.not_filtered} >> {output.not_filtered}
        """

rule find_rnap_domains:
    input:
        faa = rules.concatenate_prots.output.not_filtered,
        hmm = os.path.join(config['domains_dir'], 'pfam_rnaps.hmm')
    output:
        txt = os.path.join(config['domains_dir'], 'pfam_rnaps.txt')
    conda:
        os.path.join(config['envs'], 'hmmer.yml')
    threads: config['maxthreads']
    params: dom_e = config['hmmer_eval_thres']
    shell:
        """
        hmmsearch  --cpu {threads} \
        --noali --notextw -E {params.dom_e} --domE {params.dom_e} --tblout {output} \
        {input.hmm} {input.faa}
        """

rule convert_hmmer_txt_to_tsv:
    input:
        txt = rules.find_rnap_domains.output.txt
    output:
        tsv = os.path.join(config['domains_tables_dir'], 'pfam_rnaps.txt')
    params: script = os.path.join(config['scripts'], 'parsehmm.py')
    shell:
        """
        python {params.script} --infile {input.txt} --outfile {output.tsv}
        """


rule faa_to_gff:
    input:
        faa = rules.concatenate_prots.output.not_filtered,  # os.path.join(config['annotation_dir_prot'], '{sample}.faa'),
        domains_tsv = rules.convert_hmmer_txt_to_tsv.output.tsv
    output:
        os.path.join(config['annotation_dir'], 'concatenated.gff')
    conda:
        os.path.join(config['envs'], 'biopython.yml')
    params:
        script = os.path.join(config['scripts'], 'create_proteome_gff.py'),
        s = strandness,
        mode = 'pfam'
    shell:
        """
        python {params.script} -f {input.faa} -t {input.domains_tsv} -m {params.mode}  -o {output} {params.s}
        """


rule select_target_seqids:
    input:
        os.path.join(config['annotation_dir'], 'concatenated.gff')
    output:
        os.path.join(config['meta'], 'meta_target_contigs.txt')
    shell:
        """
        cat {input} | grep "PFAM" | cut -f 1 | sort -u > {output}
        """


rule split_target_sequences_into_separate_files:
    input:
        sequences = lambda wc: os.path.join(config['db_fastas'], f'{wc.db}', f'{paths[wc.db]}'),
        filters = os.path.join(config['meta'], 'meta_target_contigs.txt')
    output:
        directory(os.path.join(config['filtered_contigs_dir'], '{db}'))
    params:
        script=os.path.join(config['scripts'], 'split_fastas_meta.py')
    conda: os.path.join(config['envs'], 'biopython.yml')
    shell:
        """
        python {params.script} -i {input.sequences} -o {output} -f {input.filters} -s 1 -l 0
        """

rule mv_target_sequences:
    input:
        expand(os.path.join(config['filtered_contigs_dir'], '{db}'), db=databases)
    output:
        directory(config['genomes_source'])
    shell:
        """
        mkdir {output};
        array=({input})
        for dir in ${{array[*]}}
          do cp ${{dir}}/* {output}/
        done
        """