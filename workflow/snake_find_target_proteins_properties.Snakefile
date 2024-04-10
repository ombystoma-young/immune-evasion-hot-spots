import os

configfile: 'config_autographiviridae.yaml'

os.makedirs(config['target_prot_dom_db_dir'], exist_ok=True)


rule all:
    input:
        os.path.join(config['domains_tables_dir'], "target_proteins_phrogs_with_descr.tsv"),
        os.path.join(config['domains_tables_dir'],"pfam_target.tsv"),
        os.path.join(config['domains_tables_dir'], "apis_whole_with_descr.tsv"),
        os.path.join(config['target_dir'], "target_physical_char.tsv")



# FIND DOMAINS IN TARGET GENES (for all, not only representatives of clusters)
rule search_phrog:
    input:
        db = os.path.join(config['target_prot_db_dir'], 'target_proteins')
    output:
        search = os.path.join(config['target_prot_dom_db_dir'], 'target_proteins_phrogs.0')
    params:
            phrog_db = os.path.join(config['phrog_db'], 'phrogs_mmseqs_db', 'phrogs_profile_db'),
            temp_dir = os.path.join(config['mmseqs_temp'], 'tmp_phrog'),
            out_path =  os.path.join(config['target_prot_dom_db_dir'], 'target_proteins_phrogs'),
            maxmem = config['maxmem'],
            maxram = config['maxram']
    threads: config['maxthreads']
    shell:
        """
        mmseqs search \
        --split-memory-limit {params.maxram} --disk-space-limit {params.maxmem} --threads {threads} \
        -s 7 \
        {params.phrog_db} {input} {params.out_path} {params.temp_dir}
        """


# write results
rule write_results:
    input:
        search = os.path.join(config['target_prot_dom_db_dir'], 'target_proteins_phrogs.0'),
        db = os.path.join(config['target_prot_db_dir'], 'target_proteins')
    output:
        os.path.join(config['domains_tables_dir'], "target_proteins_phrogs.tsv")
    params:
            phrog_db = os.path.join(config['phrog_db'],  'phrogs_mmseqs_db', 'phrogs_profile_db'),
            input_path = os.path.join(config['target_prot_dom_db_dir'], 'target_proteins_phrogs')
    shell:
        """
        mmseqs createtsv {params.phrog_db} {input.db} {params.input_path} {output}
        """


rule add_descr_phrogs:
    input:
        os.path.join(config['domains_tables_dir'], "target_proteins_phrogs.tsv")
    output:
        os.path.join(config['domains_tables_dir'], "target_proteins_phrogs_with_descr.tsv")
    params:
        script = os.path.join(config['scripts'], 'add_domains_descr.py'),
        desc = config['phrog_desc']
    shell:
        """
        python {params.script} -d {params.desc} -i {input} -o {output} -t phrog
        """


rule search_pfam_domains:
    input:
        faa = os.path.join(config['target_dir'], "target.faa")
    output:
        os.path.join(config['domains_tables_dir'], "pfam_target.txt")
    params:
        path = config['pfam_path'],
        eval = config['hmmer_eval_thres']
    conda:
        os.path.join(config['envs'], 'hmmer.yml')
    shell:
         """    
         hmmsearch --noali --notextw -E {params.eval} --domE {params.eval} --tblout {output} {params.path} {input.faa}
         """


rule search_apis_domains:
    input:
        faa = os.path.join(config['annotation_dir'], "concatenated.faa") if config['dbapis_everywhere']
        else os.path.join(config['target_dir'], "target.faa")
    output:
        os.path.join(config['domains_tables_dir'], "apis_whole.txt")
    params:
        path = config['dbapis_path'],
        eval = config['hmmer_eval_thres']
    conda:
        os.path.join(config['envs'], 'hmmer.yml')
    shell:
         """    
         hmmsearch --noali --notextw -E {params.eval} --domE {params.eval} --tblout {output} {params.path} {input.faa}
         """


rule hmmer_res_txt2tsv:
    input:
        pfam=os.path.join(config['domains_tables_dir'], "pfam_target.txt"),
        dbapis=os.path.join(config['domains_tables_dir'], "apis_whole.txt")
    output:
        pfam=os.path.join(config['domains_tables_dir'],"pfam_target.tsv"),
        dbapis=os.path.join(config['domains_tables_dir'],"apis_whole.tsv")
    params:
        script = os.path.join(config['scripts'], 'parsehmm.py')
    shell:
        """
        python {params.script} --infile {input.pfam} --outfile {output.pfam}
        python {params.script} --infile {input.dbapis} --outfile {output.dbapis} --extractbest
        """


rule add_descr_dbapis:
    input:
        os.path.join(config['domains_tables_dir'], "apis_whole.tsv")
    output:
        os.path.join(config['domains_tables_dir'], "apis_whole_with_descr.tsv")
    params:
        script = os.path.join(config['scripts'], 'add_domains_descr.py'),
        desc = config['dbapis_desc']
    shell:
        """
        python {params.script} -d {params.desc} -i {input} -o {output} -t dbapis
        """


rule calculate_phys_characteristics:
    input:
        os.path.join(config['target_dir'], "target.faa")
    output:
        os.path.join(config['target_dir'], "target_physical_char.tsv")
    params:
        script = os.path.join(config['scripts'], 'extract_phys_chars.py')
    conda: os.path.join(config['envs'], 'biopython.yml')
    shell:
        """
        python {params.script} -i {input} -o {output}
        """
