import os

assemblies_dir = 'ncbi_dataset/data'
meta_dir = 'metadata'

blast_dir = 'blasted'
blast_db_dir = os.path.join(blast_dir, 'blastdb')
db_name = 'phages_genomes_concat'
features_dir = 'annotation'
prokka_dir = os.path.join(features_dir, 'prokka')
pharokka_dir = os.path.join(features_dir, 'pharokka')
cluster_dir = 'clusterization'
cluster_prot_dir = 'protein_clusterization_datasets'
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
os.makedirs(prokka_dir, exist_ok=True)
os.makedirs(pharokka_dir, exist_ok=True)
os.makedirs(cluster_dir, exist_ok=True)
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

genomes = os.listdir(assemblies_dir)
genomes.remove('assembly_data_report.jsonl')
genomes.remove('dataset_catalog.json')

# drop ids, soft filtering
drops = os.path.join(meta_dir, 'genomes_to_drop.tsv')
with open(drops, 'r') as drop_list:
    for line in drop_list:
        row = line.strip().split(',')
        if row[2] != '"AssemblyID"':
            id_drop = row[2][1:-1]
            genomes.remove(id_drop)

#  list of fna files for concatenation:
genome_fasta_s = []
for subdir in os.walk(assemblies_dir):
    if subdir[0] != assemblies_dir:
        genome_fasta_s.extend(["/".join([subdir[0], file]) for file in subdir[2] if file.endswith('fna')])

# remove genome ids for fasta after soft filtering
for genome_fasta in genome_fasta_s:
    genome_ = genome_fasta.split('/')[2]
    if genome_ not in genomes:
        genome_fasta_s.remove(genome_fasta)

genomes.sort()
genome_fasta_s.sort()

datasets = [f'dataset_{i}' for i in range(1, 6)]

rule all:
    input:
        os.path.join(blast_dir, 'polymerases_tblastn.tsv'),
        expand(os.path.join(prokka_dir, '{genome}', '{genome}.gff'), zip, genome=genomes),
        os.path.join(upstream_dir,'representative_genomes.gff'),
        os.path.join(upstream_dir,'upstream.gff'),
        os.path.join(domain_tables_dir, "upstream_domains.tsv"),
        expand(os.path.join(cluster_prot_dir,'upstream_proteins_{dataset}_clu.faa'), dataset=datasets),
        expand(os.path.join(results_dir, 'upstreams_with_clusters_{dataset}.gff'), dataset=datasets)
        # os.path.join(promoters_dir,'all_intergenic_with_length.tsv'),
       # os.path.join(upstream_dir,'tdr_pol_dist.tsv'),
       # os.path.join(upstream_dir,'polymerases.mafft.trim.faa'),
       # os.path.join(meta_dir,'not_all_genomic.gff'),
       # os.path.join(upstream_dir,'tree_far_','polymerases_far.iqtree.treefile')
        # os.path.join(promoters_dir, 'promoters_blastn.tsv'),   # some troubles with length of sequence
        #os.path.join(cluster_prot_dir, 'upstream_proteins_clu.faa'),
        #os.path.join(upstream_dir, 'meta_upstream.gff'),
        #os.path.join(results_dir, 'freq_repres_proteins.txt'),
        #os.path.join(results_dir,'upstreams_with_clusters.gff'),
        #os.path.join(results_dir,'upstream_proteins_clu_wide_seq_filtered.tsv')

# BLOCK PROTEIN SEQUENCES CLUSTERIZATION

rule create_prot_mmseq_db:
    input:
        os.path.join(upstream_dir, 'upstream_{dataset}.faa')
    output:
        os.path.join(cluster_prot_dir, 'upstream_proteins_{dataset}', 'upstream_proteins_{dataset}.fnaDB')
    conda:
        'envs/mmseq2.yml'
    shell:
        """
        mmseqs createdb {input} {output} --shuffle
        """

rule cluster_prot:  # лажа с такими параметрами
    input:
        db = os.path.join(cluster_prot_dir, 'upstream_proteins_{dataset}', 'upstream_proteins_{dataset}.fnaDB')
    output:
        clu = os.path.join(cluster_prot_dir, 'upstream_proteins_{dataset}', 'upstream_proteins_{dataset}_clu.0')
    conda:
        'envs/mmseq2.yml'
    params: clu = os.path.join(cluster_prot_dir, 'upstream_proteins_{dataset}', 'upstream_proteins_{dataset}_clu')
    threads: 10
    shell:
        """
        mmseqs cluster --threads {threads} --max-seqs 300 -k 6 --cluster-mode 0 \
        --alignment-mode 3 --min-seq-id 0.75 \
        --cov-mode 0 -c 0.5 \
        --split-memory-limit 7G {input} {params.clu} tmp_prot
        """

rule get_prot_clusters_tsv:
    input:
        clu = os.path.join(cluster_prot_dir, 'upstream_proteins_{dataset}', 'upstream_proteins_{dataset}_clu.0'),
        db = os.path.join(cluster_prot_dir, 'upstream_proteins_{dataset}', 'upstream_proteins_{dataset}.fnaDB')
    output:
        tsv = os.path.join(cluster_prot_dir, 'upstream_proteins_{dataset}_clu.tsv')
    params: clu=os.path.join(cluster_prot_dir, 'upstream_proteins_{dataset}', 'upstream_proteins_{dataset}_clu')
    conda:
        'envs/mmseq2.yml'
    shell:
        """
        mmseqs createtsv {input.db} {input.db} {params.clu} {output.tsv}
        """

rule get_upstream_clusters_faa:
    input:
        clu = os.path.join(cluster_prot_dir, 'upstream_proteins_{dataset}', 'upstream_proteins_{dataset}_clu.0'),
        db = os.path.join(cluster_prot_dir, 'upstream_proteins_{dataset}', 'upstream_proteins_{dataset}.fnaDB')
    output:
        faa = os.path.join(cluster_prot_dir, 'upstream_proteins_{dataset}_clu.faa')
    params: clu=os.path.join(cluster_prot_dir, 'upstream_proteins_{dataset}', 'upstream_proteins_{dataset}_clu'),
            db_path = os.path.join(cluster_prot_dir, 'upstream_proteins_{dataset}', 'DB_clu_rep')

    conda:
        'envs/mmseq2.yml'
    shell:
        """
        mmseqs createsubdb {params.clu} {input.db} {params.db_path}
        mmseqs convert2fasta {params.db_path} {output.faa}   
        """

# BLOCK ANALYSE RESULTS
rule steal:
    input:
        expand(os.path.join(assemblies_dir, '{genome}/sequence_report.jsonl'), genome=genomes)
    output:
        os.path.join(meta_dir, 'not_all_genomic.gff')
    params: dir=assemblies_dir
    shell:
        """
        cat {params.dir}/*/genomic.gff | grep -v "#" > {output}
        """


rule intersect_upstream_and_ncbi_gff:
    input:
        a = os.path.join(meta_dir, 'not_all_genomic.gff'),
        b = os.path.join(upstream_dir, 'upstream_{dataset}.gff')
    output:
        os.path.join(upstream_dir, 'meta_upstream_{dataset}.gff')
    conda: 'envs/bedtools.yml'
    shell:
        """
        bedtools intersect -s -wa -wb -a {input.a} -b {input.b} > {output}
        """


rule get_frequency_table:
    input:
        os.path.join(cluster_prot_dir, 'upstream_proteins_{dataset}_clu.tsv')
    output:
        os.path.join(results_dir, 'freq_repres_proteins_{dataset}.txt')
    shell:
        """
        cat {input} | cut -f 1 | sort | uniq -dc | sort -n > {output}
        """

rule write_freq_stat:
    input:
        os.path.join(cluster_prot_dir,'upstream_proteins_{dataset}_clu.tsv'),
        os.path.join(cluster_prot_dir,'upstream_proteins_{dataset}_clu.faa'),
        os.path.join(upstream_dir, 'representative_genomes.gff'),
        os.path.join(upstream_dir, 'meta_upstream_{dataset}.gff')
    output:
        clus_out_long = os.path.join(results_dir, 'upstream_proteins_{dataset}_clu_long.tsv'),
        clus_out_wide = os.path.join(results_dir, 'upstream_proteins_{dataset}_clu_wide.tsv'),
        clus_out_wide_seq = os.path.join(results_dir, 'upstream_proteins_{dataset}_clu_wide_seq.tsv')
    shell:
        """
        python3 scripts/steal_annotation_datasets.py {wildcards.dataset}
        """

rule sort_wide:
    input:
        clus_out_wide = os.path.join(results_dir, 'upstream_proteins_{dataset}_clu_wide.tsv')
    output:
        clus_out_wide = os.path.join(results_dir, 'upstream_proteins_{dataset}_clu_wide_sorted.tsv')
    shell:
        """
        cat {input} | sort -nrk 2 > {output}
        """


rule add_cluster_names_to_gff:
    input:
        os.path.join(results_dir, 'upstream_proteins_{dataset}_clu_wide_sorted.tsv'),
        os.path.join(results_dir, 'upstream_proteins_{dataset}_clu_long.tsv'),
        os.path.join(upstream_dir, 'upstream_{dataset}.gff')
    output:
        os.path.join(results_dir, 'upstreams_with_clusters_{dataset}.gff')
    shell:
        """
        python3 scripts/make_gff_with_clusters_datasets.py {wildcards.dataset}
        """


rule write_filtered_wide:
    input:
        gff=os.path.join(results_dir, 'upstreams_{dataset}_with_clusters.gff'),
        tsv_wide=os.path.join(results_dir, 'upstream_proteins_{dataset}_clu_wide_seq.tsv'),
        tsv_long=os.path.join(results_dir,'upstream_proteins_{dataset}_clu_long.tsv')
    output:
        pass_one = temp('genes_passes_{dataset}'),
        pass_two = temp('clusters_passes_{dataset}'),
        tsv = os.path.join(results_dir, 'upstream_proteins_{dataset}_clu_wide_seq_filtered.tsv')
    shell:
        """
        cat {input.gff} | cut -f 9 | cut -f 1 -d ";" | cut -f 2 -d "=" > {output.pass_one}
        cat {input.tsv_long} | grep -f {output.pass_one} | cut -f 1 | sort -u > {output.pass_two}
        cat {input.tsv_wide} | grep -f {output.pass_two} > {output.tsv}
        """