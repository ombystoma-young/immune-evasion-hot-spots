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
datasets_aln = [f'dataset_{i}' for i in range(1, 5)]

rule all:
    input:
        os.path.join(blast_dir, 'polymerases_tblastn.tsv'),
        expand(os.path.join(prokka_dir, '{genome}', '{genome}.gff'), zip, genome=genomes),
        os.path.join(upstream_dir,'representative_genomes.gff'),
        os.path.join(upstream_dir,'upstream_fixed.gff'),
        os.path.join(domain_tables_dir, "upstream_domains.tsv"),
        os.path.join(cluster_prot_dir,'upstream_proteins_clu.faa'),
        os.path.join(cluster_prot_dir,'upstream_proteins_clu.tsv'),
        os.path.join(results_dir, 'upstreams_with_clusters.gff'),
        os.path.join(results_dir,'upstream_proteins_clu_wide_seq_sorted.tsv'),
        os.path.join(results_dir,'upstream_proteins_clu_wide_seq_sorted_prokka.tsv')
        #os.path.join(upstream_dir, 'meta_upstream.gff'),
        #os.path.join(results_dir, 'freq_repres_proteins.txt'),
        # os.path.join(results_dir,'upstream_proteins_clu_wide_seq_filtered.tsv')

# BLOCK cut upstream regions
rule exclude_bad_tdrs:
    input:
        tsv=os.path.join(tdrs_search_dir,'TDRs_all.tsv'),
        d3=os.path.join(datasets_dir,'genomes_dataset_3.txt'),
        d4=os.path.join(datasets_dir, 'genomes_dataset_4.txt'),
        d5=os.path.join(datasets_dir,'genomes_dataset_5.txt')
    output:
        os.path.join(tdrs_search_dir, 'TDRs_modified.tsv')
    shell:
        """
        cat {input.tsv} | grep -v -f {input.d3} | grep -v -f {input.d4} | grep -v -f {input.d5} > {output}
        """

rule get_upstreams_coordinates:
    input:
        os.path.join(upstream_dir, 'representative_genomes.gff'),
        os.path.join(meta_dir, 'assembly_nuccore.tsv'),
        os.path.join(intergenic_dir, 'all_intergenic_with_length.tsv'),
        os.path.join(tdrs_search_dir,'TDRs_modified.tsv'),
        os.path.join('stats', 'genomes_gc_length.statistics'),
        os.path.join(meta_dir, 'genomes_after_curation.tsv')
    output:
        os.path.join(upstream_dir,'upstream_fixed.bed')
    script: 'scripts/getrightupstreambed.py'

rule get_upstream_genes:
    input:
        bed = os.path.join(upstream_dir, 'upstream_fixed.bed'),
        gff = os.path.join(upstream_dir, 'representative_genomes.gff')
    output:
        gff = os.path.join(upstream_dir, 'upstream_fixed.gff')
    conda: 'envs/bedtools.yml'
    shell:
        """
        cat {input.bed} | cut -f 1-6 > temp_file.bed
        bedtools intersect -a {input.gff} -b temp_file.bed -s > {output}
        rm temp_file.bed
        """

# rule get_tdr_rnap_dist:
#     input:
#         os.path.join(upstream_dir,'upstream_fixed.bed')
#     output:
#         os.path.join(upstream_dir,  'tdr_pol_dist.tsv')
#     script: 'scripts/upstream_lengths_estimator.py'


rule concat_protein_fasta:
    input:
        expand(os.path.join(prokka_dir,'{genome}','{genome}.faa'),zip,genome=genomes)
    output:
        faa = os.path.join(upstream_dir, 'all_genomes.faa')
    params: path_prokka=prokka_dir
    shell:
        """
        cat {params.path_prokka}/*/*.faa > {output} 
        """


rule get_upstream_faa:
    input:
        faa = os.path.join(upstream_dir, 'all_genomes.faa'),
        gff = os.path.join(upstream_dir, 'upstream_fixed.gff'),
    output:
        faa_total = os.path.join(upstream_dir, 'upstream.faa')
    script: 'scripts/get_upstream_proteins_faa.py'


# BLOCK PROTEIN SEQUENCES CLUSTERIZATION

rule create_prot_mmseq_db:
    input:
        os.path.join(upstream_dir, 'upstream.faa')
    output:
        os.path.join(cluster_prot_dir, 'upstream_proteins', 'upstream_proteins.fnaDB')
    conda:
        'envs/mmseq2.yml'
    shell:
        """
        mmseqs createdb {input} {output} --shuffle
        """

rule cluster_prot:
    input:
        db = os.path.join(cluster_prot_dir, 'upstream_proteins', 'upstream_proteins.fnaDB')
    output:
        clu = os.path.join(cluster_prot_dir, 'upstream_proteins_clu.0')
    conda:
        'envs/mmseq2.yml'
    params: clu = os.path.join(cluster_prot_dir, 'upstream_proteins_clu')
    threads: 10
    shell:
        """
        mmseqs cluster --threads {threads} --max-seqs 300 -k 6 --cluster-mode 1 \
        --cov-mode 0 -c 0.7 --min-seq-id 0.3 \
        --split-memory-limit 7G {input} {params.clu} tmp_prot
        """

rule get_prot_clusters_tsv:
    input:
        clu = os.path.join(cluster_prot_dir, 'upstream_proteins_clu.0'),
        db = os.path.join(cluster_prot_dir, 'upstream_proteins', 'upstream_proteins.fnaDB')
    output:
        tsv = os.path.join(cluster_prot_dir, 'upstream_proteins_clu.tsv')
    params: clu=os.path.join(cluster_prot_dir,'upstream_proteins_clu')
    conda:
        'envs/mmseq2.yml'
    shell:
        """
        mmseqs createtsv {input.db} {input.db} {params.clu} {output.tsv}
        """

rule get_upstream_clusters_faa:
    input:
        clu = os.path.join(cluster_prot_dir, 'upstream_proteins_clu.0'),
        db = os.path.join(cluster_prot_dir, 'upstream_proteins', 'upstream_proteins.fnaDB')
    output:
        faa = os.path.join(cluster_prot_dir, 'upstream_proteins_clu.faa')
    params: clu=os.path.join(cluster_prot_dir,'upstream_proteins_clu')
    conda:
        'envs/mmseq2.yml'
    shell:
        """
        mmseqs createsubdb {params.clu} {input.db} protein_clusterization/DB_clu_rep
        mmseqs convert2fasta protein_clusterization/DB_clu_rep {output.faa}   
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
        b = os.path.join(upstream_dir, 'upstream_fixed.gff')
    output:
        os.path.join(upstream_dir, 'meta_upstream.gff')
    conda: 'envs/bedtools.yml'
    shell:
        """
        bedtools intersect -s -wa -wb -a {input.a} -b {input.b} > {output}
        """


rule get_frequency_table:
    input:
        os.path.join(cluster_prot_dir, 'upstream_proteins_clu.tsv')
    output:
        os.path.join(results_dir, 'freq_repres_proteins.txt')
    shell:
        """
        cat {input} | cut -f 1 | sort | uniq -dc | sort -n > {output}
        """

rule write_freq_stat:
    input:
        os.path.join(cluster_prot_dir,'upstream_proteins_clu.tsv'),
        os.path.join(cluster_prot_dir,'upstream_proteins_clu.faa'),
        os.path.join(upstream_dir, 'representative_genomes.gff'),
        os.path.join(upstream_dir, 'meta_upstream.gff')
    output:
        clus_out_long = os.path.join(results_dir, 'upstream_proteins_clu_long.tsv'),
        clus_out_wide = os.path.join(results_dir, 'upstream_proteins_clu_wide.tsv'),
        clus_out_wide_seq = os.path.join(results_dir, 'upstream_proteins_clu_wide_seq.tsv')
    script: 'scripts/steal_annotation.py'


rule sort_wide_seq:
    input:
        clus_out_wide = os.path.join(results_dir,'upstream_proteins_clu_wide_seq.tsv')
    output:
        clus_out_wide = os.path.join(results_dir,'upstream_proteins_clu_wide_seq_sorted.tsv')
    shell:
        """
        cat {input} | sort -nrk 2 > {output}
        """

rule sort_wide:
    input:
        clus_out_wide = os.path.join(results_dir,'upstream_proteins_clu_wide.tsv')
    output:
        clus_out_wide = os.path.join(results_dir,'upstream_proteins_clu_wide_sorted.tsv')
    shell:
        """
        cat {input} | sort -nrk 2 > {output}
        """


rule only_prokka_upstreams:
    input:
        os.path.join(cluster_prot_dir,'upstream_proteins_clu.tsv'),
        os.path.join(cluster_prot_dir,'upstream_proteins_clu.faa'),
        os.path.join(upstream_dir, 'representative_genomes.gff'),
        os.path.join(upstream_dir, 'meta_upstream.gff')
    output:
        clus_out_long = os.path.join(results_dir, 'upstream_proteins_clu_long_prokka.tsv'),
        clus_out_wide = os.path.join(results_dir, 'upstream_proteins_clu_wide_prokka.tsv'),
        clus_out_wide_seq = os.path.join(results_dir, 'upstream_proteins_clu_wide_seq_prokka.tsv')
    shell:
        """
        python3 scripts/add_only_prokka_annotation.py
        """

rule sort_prokka_ann:
    input:
        os.path.join(results_dir, 'upstream_proteins_clu_wide_seq_prokka.tsv')
    output:
        os.path.join(results_dir,'upstream_proteins_clu_wide_seq_sorted_prokka.tsv')
    shell:
        """
        cat {input} | sort -rnk 2 > {output}
        """


rule filter_later_upstreams:
    input:
        os.path.join(results_dir,'upstream_proteins_clu_wide_sorted.tsv'),
        os.path.join(results_dir,'upstream_proteins_clu_long.tsv'),
        os.path.join(upstream_dir, 'upstream_fixed.gff')
    output:
        os.path.join(results_dir, 'upstreams_with_clusters.gff')
    script: 'scripts/write_gff_all_datasets.py'


rule write_filtered_wide:
    input:
        gff=os.path.join(results_dir, 'upstreams_with_clusters.gff'),
        tsv_wide=os.path.join(results_dir, 'upstream_proteins_clu_wide_seq.tsv'),
        tsv_long=os.path.join(results_dir, 'upstream_proteins_clu_long.tsv')
    output:
        pass_one = temp('genes_passes'),
        pass_two = temp('clusters_passes'),
        tsv = os.path.join(results_dir, 'upstream_proteins_clu_wide_seq_filtered.tsv')
    shell:
        """
        cat {input.gff} | cut -f 9 | cut -f 1 -d ";" | cut -f 2 -d "=" > {output.pass_one}
        cat {input.tsv_long} | grep -f {output.pass_one} | cut -f 1 | sort -u > {output.pass_two}
        cat {input.tsv_wide} | grep -f {output.pass_two} > {output.tsv}
        """


# BLOCK: DOMAINs SEARCH
rule search_domains:
    input:
        faa = os.path.join(upstream_dir, 'upstream.faa'),
        hmm = os.path.join(profiles_dir, "domains_Burstein.hmm")
    output:
        os.path.join(domain_tables_dir,"upstream_domains.tsv")
    conda: 'envs/hmmer.yml'
    shell:
        """
        hmmsearch --noali --notextw -E 0.000001 --domE 0.000001 --tblout {output} {input.hmm} {input.faa}
        """

