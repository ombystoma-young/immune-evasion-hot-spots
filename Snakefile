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
pharokka_db_dir = os.path.join('metadata', 'pharokka_db')
promoters_dir = 'promoters_search'
intergenic_regions_db = os.path.join(promoters_dir, 'ig_blast_db')

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
os.makedirs(promoters_dir, exist_ok=True)
os.makedirs(intergenic_regions_db, exist_ok=True)

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


rule all:
    input:
        os.path.join(blast_dir, 'polymerases_tblastn.ssv'),
        expand(os.path.join(prokka_dir, '{genome}', '{genome}.txt'), zip, genome=genomes),
        os.path.join(cluster_dir, 'phages_genomes_concat_clu.tsv'),
        os.path.join(upstream_dir,'representative_genomes.gff'),
        os.path.join(upstream_dir,'upstream.gff'),
        # os.path.join(cluster_prot_dir, 'upstream_proteins_clu.faa'),
        # os.path.join(cluster_prot_dir,'upstream_proteins_clu.tsv'),
        os.path.join(domain_tables_dir, "upstream_domains.tsv"),
        os.path.join(promoters_dir, 'intergenic_filtered.bed'),
        os.path.join(promoters_dir, 'promoters_blastn.tsv')

# update: report about removing useless
rule update_stat:
    input:
        stats = os.path.join(meta_dir, 'stats_for_preso')
    output:
        'update.status'
    run:
        stats = str(input.stats)
        with open(stats, 'a') as stat_file:
            stat_file.write('Number of entries after soft filter\n')
            stat_file.write(f'{len(genomes)}\n')
        open(str(output[0]), 'w').close()

rule concat_fasta:
    input:
        expand('{genome_fasts}', genome_fasts=genome_fasta_s)
    output:
        os.path.join(blast_dir, 'phages_genomes_concat.fna')
    shell:
        """
        cat {input} > {output}
        """


# # BLOCK representative genomes
# rule create_mmseq_db:
#     input:
#         os.path.join(blast_dir,'phages_genomes_concat.fna')
#     output:
#         os.path.join(cluster_dir,'phages_genomes', 'phages_genomes_concat.fnaDB')
#     conda:
#         'envs/mmseq2.yml'
#     shell:
#         """
#         mmseqs createdb {input} {output} --shuffle
#         """
#
# rule clusterization:
#     input:
#         db = os.path.join(cluster_dir,'phages_genomes','phages_genomes_concat.fnaDB')
#     output:
#         clu = os.path.join(cluster_dir,'phages_genomes_concat_clu.0')
#     conda:
#         'envs/mmseq2.yml'
#     params: clu = os.path.join(cluster_dir,'phages_genomes_concat_clu')
#     threads: 8
#     shell:
#         """
#         mmseqs cluster --threads {threads} --max-seqs 300 -k 14 --split-memory-limit 5G {input} {params.clu} tmp
#         """
#
# rule get_clusters_tsv:
#     input:
#         clu = os.path.join(cluster_dir,'phages_genomes_concat_clu.0'),
#         db = os.path.join(cluster_dir,'phages_genomes','phages_genomes_concat.fnaDB')
#     output:
#         tsv = os.path.join(cluster_dir,'phages_genomes_concat_clu.tsv')
#     params: clu=os.path.join(cluster_dir,'phages_genomes_concat_clu')
#     conda:
#         'envs/mmseq2.yml'
#     shell:
#         """
#         mmseqs createtsv {input.db} {input.db} {params.clu} {output.tsv}
#         """
#
# rule get_clusters_fasta:
#     input:
#         clu = os.path.join(cluster_dir,'phages_genomes_concat_clu.0'),
#         db = os.path.join(cluster_dir,'phages_genomes','phages_genomes_concat.fnaDB')
#     output:
#         fna = os.path.join(cluster_dir,'phages_genomes_concat_clu.fna')
#     params: clu=os.path.join(cluster_dir,'phages_genomes_concat_clu')
#     conda:
#         'envs/mmseq2.yml'
#     shell:
#         """
#         mmseqs createsubdb {params.clu} {input.db} clusterization/DB_clu_rep
#         mmseqs convert2fasta clusterization/DB_clu_rep {output.fna}
#         """
#
# # update: report about removing useless
# rule update_stat_clusters:
#     input:
#         stats = os.path.join(meta_dir, 'stats_for_preso'),
#         fna = os.path.join(cluster_dir,'phages_genomes_concat_clu.fna')
#     output:
#         'clusters.status'
#     shell:
#         """
#         echo "Number of representative genomes" >> {input.stats}
#         cat {input.fna} | grep ">" | wc -l >> {input.stats}
#         touch {output}
#         """

# BLOCK search TDRs in genomes with minimap2
rule search_TDRs:
    input:
        fna_path=os.path.join(assemblies_dir, '{genome}')
    output:
        paf=os.path.join(tdrs_search_dir, 'all', '{genome}.paf')
    params: '-X -N 50 -p 0.1 -c'
    threads: 1
    conda:
        'envs/minimap2.yml'
    shell:
        """
         minimap2 {params} {input.fna_path}/*.fna {input.fna_path}/*.fna > {output.paf}
        """

# unite pafs for analysis
rule unite_pafs:
    input:
        expand(os.path.join(tdrs_search_dir, 'all', '{genome}.paf'), genome=genomes)
    output:
        os.path.join(tdrs_search_dir, 'all_phages.paf')
    threads: 1
    shell:
        """
        cat {input} > {output}
        """

# returns tsv file with TDRs
rule find_TDRs:
    input:
        os.path.join(tdrs_search_dir, 'all_phages.paf')
    output:
        tsv = os.path.join(tdrs_search_dir, 'TDRs_all.tsv')
    threads: 1
    params: min_len=99, de=0.3
    shell:
        """ 
        Rscript scripts/find_TDRs.R {input} {output} {params} 
        """


rule replace_names_fasta:
    input:
        os.path.join(cluster_dir, 'phages_genomes_concat_clu.fna')
    output:
        os.path.join(tdrs_search_dir, 'phages_genomes_modified.fna')
    run:
        inp = str(input[0])
        print(inp)
        out = str(output[0])
        with open(out,'w') as out_f:
            with open(inp, 'r') as in_f:
                for line in in_f:
                    if line.startswith('>'):
                        name = line.strip().split()[0]
                        out_f.write(f'{name}\n')
                    else:
                        out_f.write(line)


# BLOCK blast RNAP
rule create_blast_db:
    input:
        os.path.join(blast_dir, 'phages_genomes_concat.fna')
    output:
        os.path.join(blast_db_dir, 'phages_genomes_concat.nsq')
    conda:
        'envs/blast.yml'
    shell:
        """
        makeblastdb -in {input} -dbtype nucl -out {blast_db_dir}/{db_name}
        """

rule run_tblastn:
    input:
        db = os.path.join(blast_db_dir, 'phages_genomes_concat.nsq'),
        faa = os.path.join(meta_dir, 't7p07.faa')
    output:
        os.path.join(blast_dir, 'polymerases_tblastn.ssv')
    conda:
        'envs/blast.yml'
    params:
        db_path = os.path.join(blast_db_dir, 'phages_genomes_concat'),
        fmt = '6 qseqid sseqid pident nident length mismatch gapopen qstart qend sstart send evalue bitscore',
        eval = 0.05
    threads: 8
    shell:
        """
        tblastn -query {input.faa}  -db {params.db_path}  \
        -num_threads {threads} -out {output} \
        -outfmt "{params.fmt}" \
        -evalue {params.eval}
        """

# BLOCK annotate phage sequences
rule prokka_annotation:
    input:
        os.path.join(assemblies_dir, '{genome}')
    output:
        os.path.join(prokka_dir, '{genome}/{genome}.txt')
    threads: 8
    conda: 'envs/prokka.yml'
    shell:
        """
            prokka  --force --cpus {threads} --locustag {wildcards.genome} \
            --rfam --kingdom Viruses \
            --outdir {prokka_dir}/{wildcards.genome} --prefix {wildcards.genome} \
            {input}/*.fna 
        """

rule pharokka_setup_db:
    output:
        os.path.join(pharokka_db_dir, 'CARD'),
        os.path.join(pharokka_db_dir, 'phrogs_db'),
        os.path.join(pharokka_db_dir, 'vfdb')
    conda:
        "envs/pharokka.yml"
    shell:
        """
         install_databases.py -o {pharokka_db_dir}
        """

rule pharokka_annotation:
    input:
        fna=os.path.join(assemblies_dir,'{genome}'),
        db=os.path.join(pharokka_db_dir, 'vfdb')
    output:
        os.path.join(pharokka_dir, '{genome}/{genome}.gff')
    threads: 5
    conda: 'envs/pharokka.yml'
    shell:
        """
            pharokka.py -f -p {wildcards.genome} -l {wildcards.genome} -d $PWD/{pharokka_db_dir} -t {threads} -i {input.fna}/*.fna -o {pharokka_dir}/{wildcards.genome}
        """

# BLOCK cut upstream regions
# concatenate gffs
rule concat_anno:
    input:
        expand(os.path.join(prokka_dir, '{genome}', '{genome}.gff'), zip, genome=genomes)
    output:
        os.path.join(upstream_dir, 'all_genomes.gff')
    shell:
        """
        cat {input} > {output}
        """

rule get_anno_for_representative:
    input:
        gff=os.path.join(upstream_dir, 'all_genomes.gff')
    output:
        os.path.join(upstream_dir, 'representative_genomes.gff')
    shell:
        """
        cat {input.gff} | grep "CDS" > {output}
        """

rule get_upstreams_coordinated:
    input:
        os.path.join(upstream_dir, 'representative_genomes.gff'),
        tsv= os.path.join(tdrs_search_dir,'TDRs_all.tsv')
    output:
        os.path.join(upstream_dir,'upstream.bed')
    script: 'scripts/getupstreambed.py'

rule get_upstream_genes:
    input:
        bed = os.path.join(upstream_dir,'upstream.bed'),
        gff = os.path.join(upstream_dir,'representative_genomes.gff')
    output:
        gff = os.path.join(upstream_dir, 'upstream.gff')
    conda: 'envs/bedtools.yml'
    shell:
        """
        cat {input.bed} | cut -f 1-6 > temp_file.bed
        bedtools intersect -a {input.gff} -b temp_file.bed -s > {output}
        rm temp_file.bed
        """

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
        faa = os.path.join(upstream_dir,'all_genomes.faa'),
        gff = os.path.join(upstream_dir,'upstream.gff')
    output:
        faa = os.path.join(upstream_dir,'upstream.faa')
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
    threads: 8
    shell:
        """
        mmseqs cluster --threads {threads} --max-seqs 300 -k 6 --split-memory-limit 7G {input} {params.clu} tmp_prot
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

# BLOCK: DOMAINs SEARCH

rule search_domains:
    input:
        faa = os.path.join(upstream_dir, 'upstream.faa'),
        hmm = os.path.join(profiles_dir, "domains_Burstein.hmm")
    output:
        os.path.join(domain_tables_dir,"upstream_domains.tsv")
    shell:
        """
        hmmsearch --noali --notextw -E 0.000001 --domE 0.000001 --tblout {output} {input.hmm} {input.faa}
        """

# BLOCK promoters search

rule concat_seq_reports:
    input:
        expand(os.path.join(assemblies_dir, '{genome}', 'sequence_report.jsonl'), genome=genomes)
    output:
        os.path.join(promoters_dir, 'all_sequence_report.jsonl')
    params: dir=assemblies_dir
    shell:
        """
        cat {params.dir}/*/sequence_report.jsonl > {output}
        """

rule get_nuccore_names_bed:
    input:
        os.path.join(promoters_dir, 'all_sequence_report.jsonl')
    output:
        bed = os.path.join(promoters_dir, 'chromosome_lengths.bed')
    script: 'scripts/get_chomosome_lengths.py'

rule subtract_intergenic_regions:
    input:
        bed = os.path.join(promoters_dir,'chromosome_lengths.bed'),
        gff = os.path.join(upstream_dir, 'representative_genomes.gff')
    output:
        os.path.join(promoters_dir, 'all_intergenic.bed')
    conda: 'envs/bedtools.yml'
    shell:
        """
        bedtools subtract -a {input.bed} -b {input.gff} > {output}
        """

rule filter_small_intergenic_regions:
    input:
        os.path.join(promoters_dir,'all_intergenic.bed')
    output:
        temp(os.path.join(promoters_dir,'_intergenic_filtered.bed'))
    shell:
        """
        cat {input} | awk -v FS="\t" -v OFS="\t" '$3-$2>349 {{print $0}}' > {output}
        """

rule filter_representative_intergenic_regions:
    input:
        bed=os.path.join(promoters_dir, '_intergenic_filtered.bed'),
        gff=os.path.join(upstream_dir, 'representative_genomes.gff')
    output:
        os.path.join(promoters_dir,'intergenic_filtered.bed')
    shell:
        """
        cat {input.gff} | cut -f 1 | sort -u > ids
        cat {input.bed} | grep -f ids > {output}
        rm ids
        """

rule get_fa_intergenic_regions:
    input:
        bed = os.path.join(promoters_dir, 'intergenic_filtered.bed'),
        fna = os.path.join(blast_dir, 'phages_genomes_concat.fna')
    output:
        fna = os.path.join(promoters_dir, 'intergenic.fna')
    conda: 'envs/bedtools.yml'
    shell:
        """
        bedtools getfasta -fi {input.fna} -bed {input.bed} > {output} 
        """

rule create_db_intergenic:
    input:
        fna = os.path.join(promoters_dir,'intergenic.fna')
    output:
        os.path.join(intergenic_regions_db, 'intergenic.nsq')
    params: db_path= os.path.join(intergenic_regions_db, 'intergenic')
    conda:
        'envs/blast.yml'
    shell:
        """
        makeblastdb -in {input} -dbtype nucl -out {params.db_path}
        """

rule run_blastn:
    input:
        db = os.path.join(intergenic_regions_db, 'intergenic.nsq'),
        faa = os.path.join(meta_dir, 'promoters.fna')
    output:
        os.path.join(promoters_dir, 'promoters_blastn.tsv')
    conda:
        'envs/blast.yml'
    params:
        db_path=os.path.join(intergenic_regions_db, 'intergenic'),
        fmt='6 qseqid sseqid pident nident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore',
        eval=0.005
    threads: 8
    shell:
        """
        blastn -query {input.faa}  -db {params.db_path}  \
        -task megablast \
        -num_threads {threads} -out {output} \
        -outfmt "{params.fmt}" \
        -evalue {params.eval}
        """