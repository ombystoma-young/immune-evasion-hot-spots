import os

configfile: 'config_autographiviridae.yaml'

os.makedirs(config['similarity_dir'], exist_ok=True)


rule all:
    input:
        os.path.join(config['similarity_dir'], 'res_table_short_within_communities.tsv'),
        #os.path.join(config['similarity_dir'], 'loci_communities.tsv'),
        os.path.join(config['similarity_dir'], 'target_with_clusters_phrogs_within_communities.gff')


rule filter_onegene_contigs:
    input:
        os.path.join(config['target_dir'], 'target_with_clusters.gff')
    output:
        os.path.join(config['similarity_dir'], 'target_with_clusters_no_onegene_contigs.gff')
    params:
        script = os.path.join(config['scripts'], 'filter_onegene_contigs.py')
    conda: os.path.join(config['envs'], 'num_sci_py.yml')
    shell:
        """
        python {params.script} -g {input} -o {output}
        """


rule permutation_test_2_define_min_score:
    input:
        os.path.join(config['similarity_dir'], 'target_with_clusters_no_onegene_contigs.gff')
    output:
        directory(config['permutations_dir'])
    params:
        script = os.path.join(config['scripts'], 'permutation_test.py'),
        n_permutations = config['num_permutations']
    threads: config['maxthreads']
    conda: os.path.join(config['envs'], 'num_sci_py.yml')
    shell:
        """
        python3 {params.script} -i {input} -t {threads} -n {params.n_permutations} -o {output} 
        """


rule find_edges_from_scores:
    input:
        gff = os.path.join(config['similarity_dir'], 'target_with_clusters_no_onegene_contigs.gff'),
        p = config['permutations_dir']
    output:
        os.path.join(config['similarity_dir'], 'edge_list_loci.tsv')
    params:
        script = os.path.join(config['scripts'], 'find_loci_similarity.py'),
        n_permutations = config['num_permutations']
    conda: os.path.join(config['envs'], 'num_sci_py.yml')
    shell:
        """
        python3 {params.script} -g {input.gff} -p {input.p} -n {params.n_permutations} -o {output}
        """


rule find_communities_of_loci:
    input:
        os.path.join(config['similarity_dir'], 'edge_list_loci.tsv')
    output:
        os.path.join(config['similarity_dir'], 'loci_communities.tsv')
    conda:
        os.path.join(config['envs'], 'clusterone.yml')
    params:
        path=config['cl1_path'],
        haircut_score=config['haircut_score_cl1']
    shell:
        """
        java -jar {params.path} {input} -s 2 -d 0.25 --fluff --haircut {params.haircut_score} > {output}
        """


rule find_loci_similarity_info:
    input:
        c = os.path.join(config['similarity_dir'], 'loci_communities.tsv'),
        gff = os.path.join(config['similarity_dir'], 'target_with_clusters_no_onegene_contigs.gff'),
        edges = os.path.join(config['similarity_dir'], 'edge_list_loci.tsv')
    output:
        ids = os.path.join(config['similarity_dir'], 'communities_with_known_adgs_htgs.id'),
        edges = os.path.join(config['similarity_dir'], 'edge_list_loci_updated.tsv'),
        communities = os.path.join(config['similarity_dir'], 'parsed_communities_info.tsv'),
        sizes = os.path.join(config['similarity_dir'], 'loci_sizes.tsv'),
        has_adgs = os.path.join(config['similarity_dir'], 'has_adgs.tsv')
    params:
        script = os.path.join(config['scripts'], 'extract_loci_similarity_info.py'),
        out_dir = config['similarity_dir'],
        markers = config['clans']
    conda: os.path.join(config['envs'], 'num_sci_py.yml')
    shell:
        """
        python {params.script} -c {input.c} -g {input.gff} -e {input.edges} -m {params.markers} -o {params.out_dir}
        """


# TODO: re-write this rule (now just copied from notebook)
rule aggregate_results_for_selected_communities:
    input:
        long = os.path.join(config['target_clans_concat_dir'], 'res_table_long.tsv'),
        ids = os.path.join(config['similarity_dir'], 'communities_with_known_adgs_htgs.id')
    output:
        grouped = os.path.join(config['similarity_dir'], 'res_table_short_within_communities.tsv')
    params:
        script = os.path.join(config['scripts'], 'filter_clusters_target_loci.py')
    conda: os.path.join(config['envs'], 'num_sci_py.yml')
    shell:
        """
        python {params.script} --long {input.long} --ids {input.ids} --output {output.grouped}
        """


rule filter_gff_selected_communities:
    input:
        ids = os.path.join(config['similarity_dir'], 'communities_with_known_adgs_htgs.id'),
        gff = os.path.join(config['target_dir'], 'target_with_clusters_phrogs.gff')
    output:
        gff = os.path.join(config['similarity_dir'], 'target_with_clusters_phrogs_within_communities.gff')
    shell:
        """
        cat {input.gff} | grep -f {input.ids} > {output}
        """
