import os


# define directories
source_data_dir = 'psi_blast_more_genomes'
temp_requests_dir = 'some_ids'
metadata_dir = 'metadata'
assemblies_dir = 'ncbi_dataset/data'

os.makedirs(temp_requests_dir, exist_ok=True)
os.makedirs(metadata_dir, exist_ok=True)


rule all:
    input:
        os.path.join(assemblies_dir, "assembly_data_report.jsonl")
        # os.path.join(temp_requests_dir, 'assembly_ids.txt')


# get protein ids for obtaining genome id
rule get_homolog_ids:
    input:
        os.path.join(source_data_dir, 'psi_blast_round_three_attempt_3.csv')
    output:
        temp(os.path.join(temp_requests_dir, 'protein_ids.txt'))
    shell:
        """
        cat {input} | awk -v FS="," -v OFS=","  '$4 >749 {{print $0}}'|  cut -f 2 -d "," > {output}
        """

# extract assembly IDs and species names for this phages
rule get_phage_assembly_metadata:
    input:
        os.path.join(temp_requests_dir, 'protein_ids.txt')
    output:
        os.path.join(metadata_dir, 'assembly_ids.txt')
    script: "../scripts/get_assembly_metadata.sh"

# extract assembly IDs for this phages
rule get_phage_assembly_IDs:
    input:
        os.path.join(metadata_dir, 'assembly_ids.txt')
    output:
        temp(os.path.join(temp_requests_dir, 'assembly_ids.txt'))
    shell:
        """
        cat {input} | cut -f 2 | sed '/^$/d' > {output}
        """

# use ncbi datasets to download sequence data
rule get_genomes:
    input:
        os.path.join(temp_requests_dir, 'assembly_ids.txt')
    output:
        temp("ncbi_dataset.zip")
    conda:
        "envs/datasets.yml"
    shell:
        """
        datasets download genome accession --include genome,gff3,seq-report --inputfile {input}
        """

rule unzip_data:
    input:
        "ncbi_dataset.zip"
    output:
        os.path.join(assemblies_dir, "assembly_data_report.jsonl")
        #'README.md'
    shell:
        """
        unzip {input}
        """