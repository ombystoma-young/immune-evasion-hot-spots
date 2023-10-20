import os


proteomes_dir = os.path.join('annotation', 'prokka')
alignments_dir = os.path.join('antidefence_trees', 'trimmed_alignments')
meta_dir = 'metadata'
# DEFINE HMM-PROFILES DIR
profiles_dir = 'built_hmm_profiles'

# DEFINE OUTPUT DIR
enrichment_dir = 'test_enrichment'
text_hmm_res_dir = os.path.join(enrichment_dir, 'raw_results')
hmm_res_dir = os.path.join(enrichment_dir, 'tsvs')

os.makedirs(profiles_dir, exist_ok=True)
os.makedirs(enrichment_dir, exist_ok=True)
os.makedirs(text_hmm_res_dir, exist_ok=True)
os.makedirs(hmm_res_dir, exist_ok=True)

# DEFINE SAMPLES
samples = ['ocr', 'samase']
proteomes = os.listdir(proteomes_dir)

# DEFINE PARAMETERS
hmm_e_thres = 1e-3

# WORKFLOW
rule all:
    input:
        expand([os.path.join(enrichment_dir, 'new_{sample}s_genomes.txt'),
                os.path.join(enrichment_dir, 'known_{sample}s_genomes.txt')],
            sample = samples),
        os.path.join(text_hmm_res_dir, 'additional_samase_res.txt')


rule build_hmm:
    input:
        os.path.join(alignments_dir, 'trimmed_{sample}.mafft.faa')
    output:
        os.path.join(profiles_dir, '{sample}.hmm')
    conda:
        'envs/hmmer.yml'
    threads: 10
    shell:
        """
        hmmbuild --cpu {threads} {output} {input}
        """

rule concat_proteomes:
    input:
        expand(os.path.join(proteomes_dir, '{proteome}', '{proteome}.faa'), proteome = proteomes)
    output:
        os.path.join(enrichment_dir, 'proteomes.faa')
    shell:
        """
        cat {input} > {output}
        """

rule hmm_search:
    input:
        faa = os.path.join(enrichment_dir, 'proteomes.faa'),
        hmm = os.path.join(profiles_dir, '{sample}.hmm')
    output:
        txt = os.path.join(text_hmm_res_dir, 'profile_{sample}_res.txt')
    params:
        thres = hmm_e_thres
    conda:
        'envs/hmmer.yml'
    threads: 10
    shell:
        """
        hmmsearch --cpu {threads} --noali --tblout {output.txt} -E {params.thres} --domE {params.thres} {input.hmm} {input.faa}
        """

rule convert_to_csv:
    input:
        os.path.join(text_hmm_res_dir, 'profile_{sample}_res.txt')
    output:
        os.path.join(hmm_res_dir, '{sample}s.tsv')
    shell:
        """
        python3 scripts/parsehmm.py --infile {input}  --outfile {output} 
        """

rule extract_known:
    input:
        os.path.join('antidefence_trees', 'upstreams_{sample}.gff')
    output:
        os.path.join(enrichment_dir, 'known_{sample}s.txt')
    shell:
        """
        cat {input} | cut -f 9 | cut -f 1 -d ";" | cut -f 2 -d "=" > {output}
        """

rule extract_unknown:
    input:
        f = os.path.join(enrichment_dir, 'known_{sample}s.txt'),
        tsv = os.path.join(hmm_res_dir, '{sample}s.tsv')
    output:
        os.path.join(hmm_res_dir,'new_{sample}s.tsv')
    shell:
        """
        cat {input.tsv} | grep -vf {input.f} > {output}
        """

rule extract_genomes_after_curation_ids_only:
    input:
        os.path.join('metadata', 'genomes_after_curation.tsv')
    output:
        os.path.join('metadata', 'genomes_after_curation.txt')
    shell:
        """
        cat {input} | cut -f 1 > {output}
        """

rule extract_genomes_assemby_ids_curation_only:
    input:
        nuc = os.path.join('metadata', 'genomes_after_curation.txt'),
        two_ids = os.path.join('metadata', 'assembly_nuccore.tsv'),
    output:
        os.path.join('metadata', 'genomes_after_curation_assembly.txt')
    shell:
        """
        cat  {input.two_ids} | grep -f {input.nuc} | cut -f 2 > {output}
        """

rule extract_unknown_assembly_ids:
    input:
        tsv=os.path.join(hmm_res_dir, 'new_{sample}s.tsv'),
        meta=os.path.join('metadata', 'genomes_after_curation_assembly.txt')
    output:
        new = os.path.join(enrichment_dir, 'new_{sample}s_genomes.txt'),
        known = os.path.join(enrichment_dir,'known_{sample}s_genomes.txt')
    shell:
        """
        cat {input.tsv} | grep -vf {input.meta} > {output.new}
        cat {input.tsv} | grep -f {input.meta} > {output.known} 
        """

rule test_samase_sensitivity:
    input:
        hmm = expand(os.path.join(profiles_dir, '{sample}.hmm'), sample=['samase']),
        faa = os.path.join(meta_dir, 'additional_samase.faa')
    output:
        txt=os.path.join(text_hmm_res_dir, 'additional_samase_res.txt')
    params:
        thres = hmm_e_thres
    conda:
        'envs/hmmer.yml'
    threads: 10
    shell:
        """
        hmmsearch --cpu {threads} --noali --tblout {output.txt} -E {params.thres} --domE {params.thres} {input.hmm} {input.faa}
        """
