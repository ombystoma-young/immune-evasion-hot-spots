import os

assemblies_dir = 'ncbi_dataset/data'
meta_dir = 'metadata'
taxonomy_dir = os.path.join('/', 'home', 'oxalotl', 'Tools', 'ncbi_folder')

blast_dir = 'blasted'
blast_db_dir = os.path.join(blast_dir, 'blastdb')
db_name = 'phages_genomes_concat'
cluster_prot_dir = 'protein_clusterization'
cluster_prot_by_dataset_dir = os.path.join(cluster_prot_dir, 'datasets')
pharokka_db_dir = os.path.join('metadata', 'pharokka_db')
intergenic_dir = 'promoters_search'
intergenic_regions_db = os.path.join(intergenic_dir,'ig_blast_db')
results_dir = 'results'
datasets_dir = 'define_datasets'
known_ad_dir = 'antidefence_trees'
alignments_dir = os.path.join(known_ad_dir, 'alignments')
alignments_trimmed_dir = os.path.join(known_ad_dir, 'trimmed_alignments')
trees_dir = os.path.join(known_ad_dir,  'trees')

# do not rename this (used in side-scripts):
tdrs_search_dir = 'minimap2_out'
pics_dir = 'pics'
aln_dir = 'tdr_search_aln'
upstream_dir = 'upstream_search'

profiles_dir = 'domains_hmm'
domain_tables_dir = 'domain_tables'

os.makedirs(tdrs_search_dir, exist_ok=True)
os.makedirs(os.path.join(tdrs_search_dir, 'all'), exist_ok=True)
os.makedirs(blast_dir, exist_ok=True)
os.makedirs(blast_db_dir, exist_ok=True)
os.makedirs(cluster_prot_dir, exist_ok=True)
os.makedirs(aln_dir, exist_ok=True)
os.makedirs(upstream_dir, exist_ok=True)
os.makedirs(profiles_dir, exist_ok=True)
os.makedirs(domain_tables_dir, exist_ok=True)
os.makedirs(intergenic_dir, exist_ok=True)
os.makedirs(intergenic_regions_db, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)
os.makedirs(datasets_dir, exist_ok=True)
os.makedirs(known_ad_dir, exist_ok=True)
os.makedirs(alignments_dir, exist_ok=True)
os.makedirs(trees_dir, exist_ok=True)
os.makedirs(alignments_trimmed_dir, exist_ok=True)
os.makedirs(cluster_prot_by_dataset_dir, exist_ok=True)

clusters = ['ocr', 'samase']
clu_num = {'ocr': 'cluster_num=2;', 'samase': 'cluster_num=4;|cluster_num=32;'}

rule all:
    input:
        expand(os.path.join(trees_dir, '{cluster}_bootstrap.iqtree.treefile'), cluster = clusters)


rule select_cluster_representatives:
    input:
        os.path.join(results_dir, 'upstreams_with_clusters.gff')
    output:
        os.path.join(known_ad_dir, 'upstreams_{cluster}.gff')
    params:
        keyword=lambda wildcards: clu_num[f'{wildcards.cluster}']
    shell:
        """
        cat {input} | grep -P "{params.keyword}" > {output} 
        """

rule get_protein_ids:
    input:
        os.path.join(known_ad_dir, 'upstreams_{cluster}.gff')
    output:
        os.path.join(known_ad_dir, 'protein_ids_{cluster}.tsv')
    shell:
        """
        cat {input} | tr -s "\t" ";" | cut -f 1,9 -d ";" | tr -s "=" ";" | cut -f 1,3 -d ";"  | tr -s ";" "\t" > {output}
        """


# BLOCK ALIGN RNAPS FOR BEST DATASET: MAFFT + TRIMAL + IQTREE
rule get_cluster_faa:
    input:
        faa = os.path.join(upstream_dir, 'all_genomes.faa'),
        tsv = os.path.join(known_ad_dir, 'protein_ids_{cluster}.tsv'),
    output:
        faa = os.path.join(known_ad_dir,  'upsteam_{cluster}.faa')
    shell:
        """
       python scripts/get_cluster_faa.py {input.faa} {input.tsv} {output}
       """

rule align_proteins:
    input:
        os.path.join(known_ad_dir,  'upsteam_{cluster}.faa')
    output:
        os.path.join(alignments_dir,  'mafft_{cluster}.faa')
    conda: 'envs/mafft.yml'
    threads: 10
    shell:
        """
         mafft --thread {threads} --maxiterate 1000 --globalpair {input} > {output} 
        """


rule filter_alignment:
    input:
        os.path.join(alignments_dir,  'mafft_{cluster}.faa')
    output:
        os.path.join(alignments_trimmed_dir, 'trimmed_{cluster}.mafft.faa')
    conda: 'envs/trimal.yml'
    shell:
        """
        trimal -in {input} -out {output} -automated1
        """


rule build_tree_iqtree:
    input:
        os.path.join(alignments_trimmed_dir, 'trimmed_{cluster}.mafft.faa')
    output:
        os.path.join(trees_dir, '{cluster}_bootstrap.iqtree.treefile')
    params:
        pref=os.path.join(trees_dir, '{cluster}_bootstrap.iqtree'),
        model='LG+I+R6',
        ubootstrap=1000
    threads: 10
    conda: 'envs/iqtree2.yml'
    shell:
        """
        iqtree2 -nt {threads} -m {params.model} -s {input[0]} --prefix {params.pref} -B {params.ubootstrap}
        """
