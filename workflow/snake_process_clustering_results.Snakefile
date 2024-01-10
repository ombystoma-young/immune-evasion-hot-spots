import os

configfile: 'config_autographiviridae_refseq.yaml'


rule all:
    input:
        os.path.join(config['early_clans_concat_dir'], 'early_clans_edgelist.tsv'),
        os.path.join(config['early_clans_concat_dir'], 'clu2clans_parsed.tsv'),
        os.path.join(config['early_clans_concat_dir'], 'clans_wide_for_sankey.tsv'),
        os.path.join(config['early_clans_concat_dir'], "clusters_description.tsv"),
        os.path.join(config['early_clans_concat_dir'], "res_table_long.tsv"),
        os.path.join(config['early_clans_concat_dir'], "res_table.tsv"),
        os.path.join(config['upstreams_dir'],'early_with_clusters_phrogs.gff')

# EXTRACT CLANS INFO
## parse similarity
rule reformat_hhsearch_res:
    input:
        os.path.join(config['early_clans_concat_dir'], 'early_clans_info.tsv')
    output:
        os.path.join(config['early_clans_concat_dir'], 'early_clans_narrow.tsv')
    shell:
        """
        cat {input}  | sed 's/cl-//' | sed 's/|Representative=//' | cut -f 1,2,3,11 > {output}
        """

rule create_edgelist:
    input:
        os.path.join(config['early_clans_concat_dir'], 'early_clans_narrow.tsv')
    output:
        os.path.join(config['early_clans_concat_dir'], 'early_clans_edgelist.tsv')
    params:
        script = os.path.join(config['scripts'], 'create_edgelist_protein_clusters.py'),
        eval_thres = config['hhsearch_eval_thres']
    shell:
        """
        python {params.script} -i {input} -o {output} -e {params.eval_thres}
        """

## cluster clusters to clans
rule cluster_clusters_to_clans:
    input:
        os.path.join(config['early_clans_concat_dir'], 'early_clans_edgelist.tsv')
    output:
        os.path.join(config['early_clans_concat_dir'], 'clans.tsv')
    conda:
        os.path.join(config['envs'], 'clusterone.yml')
    params:
        path = config['cl1_path'],
        haircut_score = config['haircut_score_cl1']

    shell:
        """
        java -jar {params.path} {input} --fluff --haircut {params.haircut_score} > {output}
        """

rule parse_cl1_res:
    input:
        i=os.path.join(config['early_clans_concat_dir'], 'clans.tsv'),
        d=os.path.join(config['early_clans_concat_dir'], 'clusters_description.tsv'),
        e=os.path.join(config['early_clans_concat_dir'], 'early_clans_edgelist.tsv'),
    output:
        w=os.path.join(config['early_clans_concat_dir'], 'clu2clans_parsed.tsv'),
        l=os.path.join(config['early_clans_concat_dir'], 'clans_wide_for_sankey.tsv')
    params:
        script = os.path.join(config['scripts'], 'reformat_clu1_res.py')
    shell:
        """
        python {params.script} -i {input.i} -d {input.d} -e {input.e} -o {output.l} -w {output.w} 
        """

rule extract_clusters_info:
    input:
        os.path.join(config['early_clu_msa_dir'], 'msa.faDB')
    output:
        os.path.join(config['early_clans_concat_dir'], "clusters_description.tsv")
    params:
        script = os.path.join(config['scripts'], 'write_clusters_info.py'),
        in_dir = config['early_msa_unpacked_dir']
    shell:
        """
        python {params.script} -i {params.in_dir} -o {output}
        """

rule add_proteins_info:
    input:
        clu_descr = os.path.join(config['early_clans_concat_dir'], "clusters_description.tsv"),
        phrogs = os.path.join(config['domains_tables_dir'], "early_proteins_phrogs_with_descr.tsv"),
        pfam = os.path.join(config['domains_tables_dir'], "pfam_early.tsv"),
        apis = os.path.join(config['domains_tables_dir'], "apis_whole_with_descr.tsv"),
        phys_char = os.path.join(config['upstreams_dir'], "early_physical_char.tsv"),
        clans = os.path.join(config['early_clans_concat_dir'],'clu2clans_parsed.tsv')
    output:
        res_tab_long = os.path.join(config['early_clans_concat_dir'], "res_table_long.tsv"),
        res_tab = os.path.join(config['early_clans_concat_dir'], "res_table.tsv")
    params:
        script = os.path.join(config['scripts'], 'unite_protein_info_to_clusters.py')
    shell:
        """
        python {params.script} --clu {input.clu_descr} --phrogs {input.phrogs} --pfam {input.pfam} \
        --apis {input.apis} --phys {input.phys_char} --clans {input.clans} \
        --output {output.res_tab} --long {output.res_tab_long}
        """

rule select_cols_for_gff_formation:
    input:
        os.path.join(config['early_clans_concat_dir'], "res_table_long.tsv")
    output:
        os.path.join(config['early_clans_concat_dir'], "prot2famclan.tsv")
    shell:
        """
        cat {input} | awk -v FS="\t" -v OFS="\t" '{{print $2,$3,$1,$16}}' > {output} 
        """

rule add_clusters_to_gff:
    input:
        tsv = os.path.join(config['early_clans_concat_dir'], "prot2famclan.tsv"),
        gff =  os.path.join(config['upstreams_dir'], 'early.gff')
    output:
        gff = os.path.join(config['upstreams_dir'], 'early_with_clusters.gff')
    params:
        script = os.path.join(config['scripts'], 'add_clusters_info_to_gff.py')
    shell:
        """
        python {params.script} -i {input.gff} -c {input.tsv} -o {output.gff} 
        """

rule add_phrogs_descr_to_gff:
    input:
        tsv = config['phrog_desc'],
        gff = os.path.join(config['upstreams_dir'], 'early_with_clusters.gff')
    output:
        gff = os.path.join(config['upstreams_dir'], 'early_with_clusters_phrogs.gff')
    params:
        script = os.path.join(config['scripts'], 'add_phrog_description_to_gff.py')
    shell:
        """
        python {params.script} -i {input.gff} -p {input.tsv} -o {output.gff} 
        """