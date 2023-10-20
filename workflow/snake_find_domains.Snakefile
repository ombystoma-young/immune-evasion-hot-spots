import os

assemblies_dir = 'ncbi_dataset/data'
meta_dir = 'metadata'
databases_dir = os.path.join('/', 'home', 'oxalotl', 'Tools', 'databases')
pfam_dir = os.path.join(databases_dir, 'pfam_21_11_15')
phrog_dir = os.path.join(databases_dir, 'phrogs_v_4')

status_dir = 'statuses'
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
phrog_search_dir = 'phrog_search_mmseq'
phrog_seq_db_dir = 'phrog_db_mmseq'

# do not rename this (used in side-scripts):
tdrs_search_dir = 'minimap2_out'
pics_dir = 'pics'
aln_dir = 'tdr_search_aln'
upstream_dir = 'upstream_search'

profiles_dir = 'domains_hmm'
domain_tables_dir = 'domain_tables'
databases = ['phrog']
paths = {'pfam': os.path.join(pfam_dir, 'Pfam-A.hmm'),
         'phrog': os.path.join(phrog_dir, 'HMM_phrog.hmm')}


os.makedirs(status_dir, exist_ok=True)
os.makedirs(phrog_search_dir, exist_ok=True)
os.makedirs(phrog_seq_db_dir, exist_ok=True)

rule all:
    input:
        expand(os.path.join(domain_tables_dir, "{database}.tsv"), database=databases),
        os.path.join(domain_tables_dir, "phrog.tsv")

rule prepare_dbs:
    input:
        pfam = os.path.join(pfam_dir, 'Pfam-A.hmm.gz')
    output:
        pfam = os.path.join(pfam_dir, 'Pfam-A.hmm'),
        status = expand(os.path.join(status_dir, '{database}.status'), database=databases)
    params: concat_path_phrog=os.path.join(phrog_dir,'HMM_phrog')
    shell:
        """
        gunzip -k {input.pfam}  
        touch {output.status}
        """

rule hmmpress:
    input:
        status = os.path.join(status_dir, '{database}.status')
    output:
        status = os.path.join(status_dir, '{database}_hmmpress.status')
    conda:
        'envs/hmmer.yml'
    params:
        hmm = lambda wildcards: paths[f'{wildcards.database}']
    shell:
        """
        hmmpress {params.hmm}
        touch {output.status}
        """


rule search_domains:
    input:
        os.path.join(status_dir, '{database}.status'),
        faa = os.path.join(upstream_dir, "upstream.faa")
    output:
        os.path.join(domain_tables_dir, "{database}.tsv")
    params:
        path = lambda wildcards: paths[f'{wildcards.database}']
    conda:
        'envs/hmmer.yml'
    shell:
         """    
         hmmsearch --noali --notextw -E 0.0001 --domE 0.0001 --tblout {output} {params.path} {input.faa}
         """

### DONE ON SERVER:
# rule create_db_phrog:
#     input:
#         faa = os.path.join(upstream_dir, "upstream.faa")
#     output:
#         os.path.join(phrog_seq_db_dir, 'upstream')
#     conda:
#         'envs/mmseq2.yml'
#     shell:
#         """
#         mmseqs createdb {input} {output}
#         """
#
# rule search_phrog:
#     input:
#         db = os.path.join(phrog_seq_db_dir, 'upstream')
#     output:
#         search = os.path.join(phrog_search_dir, 'upstream_phrog.0')
#     conda:
#         'envs/mmseq2.yml'
#     params:
#             out = os.path.join(phrog_search_dir, 'upstream_phrog'),
#             phrog_db = os.path.join(phrog_dir, 'phrogs_mmseqs_db', 'phrogs_profile_db')
#     threads: 10
#     shell:
#         """
#         mmseqs search {params.phrog_db} {input} {params.out} tmp_phrog -s 7
#         """
#
# rule write_results:
#     input:
#         search = os.path.join(phrog_search_dir, 'upstream_phrog.0'),
#         db= os.path.join(phrog_seq_db_dir,'upstream')
#     output:
#         os.path.join(domain_tables_dir, "phrog.tsv")
#     params:
#             search = os.path.join(phrog_search_dir, 'upstream_phrog'),
#             phrog_db = os.path.join(phrog_dir,  'phrogs_mmseqs_db', 'phrogs_profile_db')
#     conda:
#         'envs/mmseq2.yml'
#     shell:
#         """
#         mmseqs createtsv {params.phrog_db} {input.db} {params.search} {output}
#         """
        
