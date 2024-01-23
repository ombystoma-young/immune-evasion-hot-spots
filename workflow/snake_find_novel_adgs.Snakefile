import os

configfile: 'config_autographiviridae_refseq.yaml'

os.makedirs(config['similarity_dir'], exist_ok=True)
os.makedirs(config['permutations_dir'], exist_ok=True)


rule all:
    input:
        os.path.join(config['similarity_dir'], 'edge_list_loci.tsv')

rule permutation_test_2_define_min_score:  # TODO: argparse
    input:
        os.path.join(config['upstreams_dir'], 'early_with_clusters.gff')
    output:
        directory(config['permutations_dir'])
    params:
        script = os.path.join(config['scripts'], 'permutation_test.py'),
        n_permutations = config['num_permutations']
    threads: 8
    shell:
        """
        python3 {params.script}
        """

rule find_edges_from_scores:  # TODO: argparse
    input:
        gff = os.path.join(config['upstreams_dir'],'early_with_clusters.gff'),
        p = expand(os.path.join(config['permutations_dir'], '{i}.pickle'),
                                i=range(config['num_permutations']))
    output:
        os.path.join(config['similarity_dir'], 'edge_list_loci.tsv')
    params:
        script = os.path.join(config['scripts'], 'find_loci_similarity.py'),
        in_dir = config['permutations_dir']
    shell:
        """
        python3 {params.script}
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
        java -jar {params.path} {input} --fluff --haircut {params.haircut_score} > {output}
        """