import os
import pandas as pd


configfile: 'config_autographiviridae.yaml'


os.makedirs(config['lineage_dir'], exist_ok=True)


rule all:
    input:
        expand(os.path.join(config['lineage_dir'], 'autographiviridae_{source}_parsed_joined_lineage.tsv'),
            source = ['phage', 'host'])


rule extract_taxonid:
    output:
        os.path.join(config['lineage_dir'], 'autographiviridae_phage.tax.ids')
    params:
        descr_path = config['ncbi_metadata'],
        script = os.path.join(config['scripts'], 'extract_taxids.py')
    shell:
        """
        python {params.script} -i {params.descr_path} -o {output}
        """


rule extract_host_name:
    output:
        os.path.join(config['lineage_dir'], 'autographiviridae_hosts.tsv')
    params:
        descr_path = config['ncbi_metadata'],
        script = os.path.join(config['scripts'], 'extract_hosts.py')
    shell:
        """
        python {params.script} -i {params.descr_path} -o {output}
        """


rule extract_host_lineage:
    input:
        os.path.join(config['lineage_dir'], 'autographiviridae_hosts.tsv')
    output:
        lin = os.path.join(config['lineage_dir'], 'autographiviridae_host.lineage'),
        tsv = os.path.join(config['lineage_dir'],'autographiviridae_host.tax.ids')
    conda:
        os.path.join(config['envs'], 'taxonkit.yml')
    params:
        db_path = config['taxdump_path']
    shell:
        """
        cat {input} | cut -f 2 | taxonkit name2taxid --data-dir {params.db_path} - > {output.tsv}
        cat {output.tsv} | cut -f 2 | taxonkit lineage - --data-dir {params.db_path} -R > {output.lin}
        """


rule extract_virus_lineage:
    input:
        os.path.join(config['lineage_dir'], 'autographiviridae_phage.tax.ids')
    output:
        os.path.join(config['lineage_dir'], 'autographiviridae_phage.lineage')
    conda:
        os.path.join(config['envs'], 'taxonkit.yml')
    params:
        db_path = config['taxdump_path']
    shell:
        """
        cat {input} | cut -f 2 | taxonkit lineage - --data-dir {params.db_path} -R > {output}
        """


rule parse_taxkit:
    input:
        os.path.join(config['lineage_dir'],  'autographiviridae_{source}.lineage')
    output:
        os.path.join(config['lineage_dir'], 'autographiviridae_{source}_parsed.lineage')
    params:
        script = os.path.join(config['scripts'], 'parse_taxkit_lineage.py')
    shell:
        """
        python3 {params.script} -i {input} -o {output} 
        """


rule join_with_ids:
    input:
        lin = os.path.join(config['lineage_dir'], 'autographiviridae_{source}_parsed.lineage'),
        ids = os.path.join(config['lineage_dir'],'autographiviridae_{source}.tax.ids')
    output:
        os.path.join(config['lineage_dir'], 'autographiviridae_{source}_parsed_joined_lineage.tsv')
    params:
        script = os.path.join(config['scripts'], 'join_tax.py'),
        hosts = lambda wc: f"-s {os.path.join(config['lineage_dir'], 'autographiviridae_hosts.tsv')}" if wc.source == 'host' else ''
    shell:
        """
        python {params.script} -i {input.lin} -t {input.ids} -o {output} {params.hosts}
        """