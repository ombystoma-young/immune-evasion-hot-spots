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
alignments_dir = os.path.join(datasets_dir, 'alignments')
trees_dir = os.path.join(datasets_dir, 'trees')

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
os.makedirs(alignments_dir, exist_ok=True)
os.makedirs(trees_dir, exist_ok=True)
os.makedirs(cluster_prot_by_dataset_dir, exist_ok=True)


rule all:
    input:
        os.path.join(trees_dir, 'polymerases_all.iqtree.treefile'),
        os.path.join(meta_dir,'genome_lineage.tsv')


# BLOCK ALIGN RNAPS FOR BEST DATASET: MAFFT + TRIMAL + IQTREE
rule get_rnap_faa:
    input:
        faa = os.path.join(upstream_dir, 'all_genomes.faa'),
        gff = os.path.join(upstream_dir, 'representative_genomes.gff'),
        tsv = os.path.join(datasets_dir,  'joined.tsv'),
        list_ = os.path.join(meta_dir, 'genomes_after_curation.tsv')
    output:
        faa = os.path.join(alignments_dir,  'polymerases_all.faa')
    shell:
        """
       python scripts/get_RNAP_faa.py all
       """

rule align_rnaps:
    input:
        os.path.join(alignments_dir,  'polymerases_all.faa')
    output:
        os.path.join(upstream_dir,'polymerases_all.mafft.faa')
    conda: 'envs/mafft.yml'
    threads: 10
    shell:
        """
         mafft --thread {threads} --maxiterate 1000 --globalpair {input} > {output} 
        """


rule filter_alignment:
    input:
        os.path.join(upstream_dir,'polymerases_all.mafft.faa')
    output:
        os.path.join(upstream_dir,'polymerases_all.mafft.trim.faa')
    conda: 'envs/trimal.yml'
    shell:
        """
        trimal -in {input} -out {output} -automated1
        """


rule build_tree_iqtree:
    input:
        os.path.join(upstream_dir,'polymerases_all.mafft.trim.faa')
    output:
        os.path.join(trees_dir, 'polymerases_all.iqtree.treefile')
    params:
        pref=os.path.join(trees_dir, 'polymerases_all.iqtree'),
        model='LG+I+R6',
        bootstrap=100
    threads: 10
    conda: 'envs/iqtree2.yml'
    shell:
        """
        iqtree2 -nt {threads} -m {params.model} -s {input[0]} --prefix {params.pref} # -b {params.bootstrap}
        """


rule get_taxons:
    input:
        genomes = os.path.join(meta_dir, 'genomes_after_curation.tsv'),
        tax = os.path.join(taxonomy_dir, 'nucl_gb.accession2taxid.gz')
    output:
        os.path.join(meta_dir, 'genome_id2taxid.tsv')
    shell:
        """
        cat {input.genomes} | cut -f 1 > temp_genomes
        zcat {input.tax} | grep -f temp_genomes > {output}
        rm temp_genomes
        """

rule get_lineage:
    input:
        os.path.join(meta_dir,'genome_id2taxid.tsv')
    output:
        os.path.join(meta_dir,'genome_lineage.tsv')
    conda:
        "envs/taxonkit.yml"  # check manual firstly, uses DOWNLOADED NCBI taxonomy db
    shell:
        """
        cat {input} | cut -f 3 | taxonkit lineage | tr -s "\t" ";" > {output}
        """
        
rule get_host_stuff:
    input:
        os.path.join(meta_dir, 'host_ids.tsv')
    output:
        os.path.join(meta_dir, 'host2taxid.tsv'),
        lin=os.path.join(meta_dir, 'host_lineage.tsv')
    conda:
        "envs/taxonkit.yml"  # check manual firstly, uses DOWNLOADED NCBI taxonomy db
    shell:
        """
        cat {input} | cut -f 2 | grep -v "NA" | sort -u > host_genus_temp.txt
        taxonkit name2taxid host_genus_temp.txt > host2taxid.tsv
        cat host2taxid.tsv | cut -f 2 | taxonkit lineage | tr -s "\t" ";" > {output.lin}
        """

