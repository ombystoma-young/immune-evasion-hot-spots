import os

configfile: 'config_autographiviridae_refseq.yaml'

os.makedirs(config['similarity_dir'], exist_ok=True)


rule all:
    input:
        os.path.join(config['similarity_dir'], 'loci_communities.tsv')
        # config['permutations_dir']


rule filter_onegene_contigs:
    input:
        os.path.join(config['upstreams_dir'], 'early_with_clusters.gff')
    output:
        os.path.join(config['similarity_dir'], 'early_with_clusters_no_onegene_contigs.gff')
    params:
        script = os.path.join(config['scripts'], 'filter_onegene_contigs.py')
    conda: os.path.join(config['envs'], 'num_sci_py.yml')
    shell:
        """
        python {params.script} -g {input} -o {output}
        """


rule permutation_test_2_define_min_score:
    input:
        os.path.join(config['similarity_dir'], 'early_with_clusters_no_onegene_contigs.gff')
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
        gff = os.path.join(config['upstreams_dir'], 'early_with_clusters_no_onegene_contigs.gff'),
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