import os

configfile: 'config_autographiviridae.yaml'

os.makedirs(config['intergenic_dir'], exist_ok=True)
os.makedirs(config['target_dir'], exist_ok=True)
os.makedirs(config['pics'], exist_ok=True)

# DEFINE PARAMETERS
if config['genomes_type'] == 'refseq':
    meta_mode = ''
    os.makedirs(config['repeats_dir'], exist_ok=True)
    os.makedirs(config['repeats_per_genome_dir'], exist_ok=True)
else:
    meta_mode = '--meta'

if config['filter_strand']:
    strandness = '-s'
else:
    strandness = ''

# READ SAMPLES
samples = []
for sample in os.listdir(config['genomes_source']):
    if config['genomes_type'] == 'refseq':
        if sample.startswith('GC'):
            samples.append(sample)
    else:
        if not sample.startswith('concat') and not sample.startswith('.snakemake'):
            samples.append(sample.replace('.fna', ''))


rule all:
    input:
        os.path.join(config['target_dir'], 'target.faa')


rule extract_RNAP_coordinates:
    input:
        gff = os.path.join(config['annotation_dir'], 'concatenated.gff'),
        meta_domains = os.path.join(config['meta'], 'rnap_phrogs.txt') if config['genomes_type'] == 'refseq' else os.path.join(config['meta'], 'rnaps_hmm_list.txt')
    output:
        os.path.join(config['annotation_dir'], 'concatenated_rnaps_only.gff')
    params: script = os.path.join(config['scripts'], 'extract_earliest_rnap.py')
    shell:
        """
        cat {input.gff} | grep -f {input.meta_domains} > temp
        python {params.script} temp {output}
        rm temp
        """


# search TDRs in genomes with minimap2
rule search_repeats:
    input:
        jsonl = os.path.join(config['genomes_source'], '{sample}', 'sequence_report.jsonl') if
            config['genomes_type'] == 'refseq' else
                os.path.join(config['genomes_source'], '{sample}.fna')
    output:
        paf = os.path.join(config['repeats_per_genome_dir'], '{sample}.paf')
    params:
        params='-X -N 50 -p 0.1 -c',
        path = os.path.join(config['genomes_source'], '{sample}', '*.fna') if
            config['genomes_type'] == 'refseq' else
                os.path.join(config['genomes_source'], '{sample}.fna')
    threads: 1
    conda:
        os.path.join(config['envs'], 'minimap2.yml')
    shell:
        """
        minimap2 {params.params} {params.path} {params.path} > {output.paf}
        """

# unite pafs for analysis
rule unite_pafs:
    input:
        expand(os.path.join(config['repeats_per_genome_dir'], '{sample}.paf'), sample=samples)
    output:
        os.path.join(config['repeats_dir'], 'tdrs.paf')
    threads: 1
    shell:
        """
        cat {input} > {output}
        """

# returns tsv file with TDRs
rule find_TDRs:
    input:
        paf = os.path.join(config['repeats_dir'], 'tdrs.paf')
    output:
        tsv = os.path.join(config['repeats_dir'], 'tdrs.tsv'),
        png = os.path.join(config['pics'], 'hist_length_TDRs.png')
    threads: 1
    params:
            min_len = config['min_tdr_len'],
            de = config['min_tdr_de'],
            script = os.path.join(config['scripts'], 'find_TDRs.R')
    conda:
        os.path.join(config['envs'], 'r-env.yml')
    shell:
        """
        Rscript {params.script} {input} {output.tsv} {output.png} {params.min_len} {params.de} 
        """

rule create_masked_low_score:
    input:
        os.path.join(config['annotation_dir'], 'concatenated.gff')
    output:
        os.path.join(config['annotation_dir'], 'concatenated_masked.gff')
    params:
        threshold = config['masking_threshold']
    shell:
        """
        cat {input} | awk -v t={params.threshold} -v FS="\t" -v OFS="\t" '$6<t {{print $0}}' > {output}
        """

rule concat_fasta:
    input:
        expand(os.path.join(config['genomes_source'], '{sample}', 'sequence_report.jsonl'), sample=samples) if
            config['genomes_type'] == 'refseq' else
                expand(os.path.join(config['genomes_source'], '{sample}.fna'), sample=samples)
    output:
        os.path.join(config['genomes_source'], 'concatenated_genomes.fna')
    params:
        paths = os.path.join(config['genomes_source'], '*', '*.fna') if
            config['genomes_type'] == 'refseq' else
                os.path.join(config['genomes_source'], '*.fna')
    shell:
        """
        cat {params.paths} > {output}
        """

rule get_lengths_bed:
    input:
        os.path.join(config['genomes_source'], 'concatenated_genomes.fna')
    output:
        bed = os.path.join(config['intergenic_dir'], 'chromosome_lengths.bed')
    params:
        script = os.path.join(config['scripts'], 'get_chomosome_lengths.py')
    conda:
        os.path.join(config['envs'], 'biopython.yml')
    shell:
        """
        python {params.script} -f {input} -b {output}
        """

rule substract_intergenic_regions:
    input:
        bed=os.path.join(config['intergenic_dir'], 'chromosome_lengths.bed'),
        gff=os.path.join(config['annotation_dir'], 'concatenated_masked.gff')
    output:
        os.path.join(config['intergenic_dir'], 'all_intergenic.bed')
    conda: 'envs/bedtools.yml'
    shell:
        """
        bedtools subtract -a {input.bed} -b {input.gff} > {output}
        """

rule get_intergenic_length:
    input:
        bed=os.path.join(config['intergenic_dir'], 'all_intergenic.bed'),
        lens = os.path.join(config['intergenic_dir'], 'chromosome_lengths.bed')
    output:
        os.path.join(config['intergenic_dir'], 'all_intergenic_with_length.tsv')
    params:
        script = os.path.join(config['scripts'], 'extract_intergenic_ends.py')
    shell:
        """
        python {params.script} -i {input.bed} -l {input.lens} -o {output}
        """


rule plot_tdr_intergenic_statistics:
    input:
        intergenics = os.path.join(config['intergenic_dir'], 'all_intergenic_with_length.tsv'),
        rnaps = os.path.join(config['annotation_dir'], 'concatenated_rnaps_only.gff'),
        tdrs = os.path.join(config['repeats_dir'],'tdrs.tsv'),
        lens = os.path.join(config['intergenic_dir'],'chromosome_lengths.bed')
    output:
        intergenic = os.path.join(config['pics'], 'diagnostic_plot_intergenic_dist_RNAP.pdf'),
        tdr = os.path.join(config['pics'], 'diagnostic_plot_tdr_RNAP.pdf'),
        all_intergenics = os.path.join(config['pics'], 'diagnostic_plot_intergenic_RNAP.pdf')
    params:
        script = os.path.join(config['scripts'], 'plot_tdr_intergenic_diagnostics.R')
    conda:
        os.path.join(config['envs'], 'r-env.yml')
    shell:
        """
        Rscript {params.script} {input.tdrs} {input.intergenics} {input.rnaps} {input.lens} \
        {output.tdr} {output.intergenic} {output.all_intergenics}
        """


rule filter_based_on_distance:
    input:
        rnaps = os.path.join(config['annotation_dir'],'concatenated_rnaps_only.gff'),
        tdrs = os.path.join(config['repeats_dir'],'tdrs.tsv'),
        intergenics = os.path.join(config['intergenic_dir'],'all_intergenic_with_length.tsv'),
        lengths = os.path.join(config['intergenic_dir'],'chromosome_lengths.bed'),
    output:
        txt = os.path.join(config['meta'], 'filtered_upstreams_nuccore.id'),
        tdrs = os.path.join(config['repeats_dir'], 'best_tdrs.tsv'),
        igs = os.path.join(config['intergenic_dir'], 'best_intergenics.tsv'),
    params:
        script = os.path.join(config['scripts'], 'split_into_datasets.R'),
        maxdisttdr = config['max_dist_to_pol_tdr'],
        maxdistinter = config['max_dist_to_pol_inter'],
        mindistinter = config['min_dist_to_pol_inter'],
        maxgen = config['max_genome_len'],
        minleninter = config['min_intergenic_len'],
        maxleninter = config['max_intergenic_len']
    conda:
        os.path.join(config['envs'], 'r-env.yml')
    shell:
        """
        Rscript {params.script} {input.tdrs} {input.intergenics} {input.rnaps} {input.lengths} {params.maxdisttdr} {params.maxdistinter} {params.mindistinter} {params.maxgen} {params.minleninter} {params.maxleninter} {output.txt} {output.tdrs} {output.igs}
        """


rule find_upstreams_coordinates:
    input:
        rnaps = os.path.join(config['annotation_dir'], 'concatenated_rnaps_only.gff'),
        tdrs = os.path.join(config['repeats_dir'], 'best_tdrs.tsv'),
        intergenics = os.path.join(config['intergenic_dir'], 'best_intergenics.tsv'),
        lengths = os.path.join(config['intergenic_dir'], 'chromosome_lengths.bed'),
        filtered = os.path.join(config['meta'], 'filtered_upstreams_nuccore.id') if not config['manual_choice_genomes'] else config['manual_choice_path']
    output:
        os.path.join(config['target_dir'], 'target.bed')
    params:
        script = os.path.join(config['scripts'], 'get_upstream_bed.py'),
        distance = config['max_dist_from_pol'],
    shell:
        """
        python {params.script} -r {input.rnaps} -t {input.tdrs} -i {input.intergenics} \
        -l {input.lengths} -f {input.filtered} -b {params.distance} -o {output}
        """

rule filter_intergenics_meta:
    input:
        rnaps = os.path.join(config['annotation_dir'], 'concatenated_rnaps_only.gff'),
        intergenics = os.path.join(config['intergenic_dir'], 'all_intergenic_with_length.tsv'),
        lengths = os.path.join(config['intergenic_dir'], 'chromosome_lengths.bed'),
    output:
        igs = os.path.join(config['intergenic_dir'], 'best_intergenics_meta.tsv')
    conda:
        os.path.join(config['envs'], 'r-env.yml')
    params:
        script = os.path.join(config['scripts'], 'split_into_datasets_meta.R'),
        maxdistinter = config['max_dist_to_pol_inter'],
        mindistinter = config['min_dist_to_pol_inter'],
        maxgen = config['max_genome_len'],
        minleninter = config['min_intergenic_len'],
        maxleninter = config['max_intergenic_len'],
        mincirclen = config['min_contig_len_circ']
    shell:
        """
        Rscript {params.script} {input.intergenics} {input.rnaps} {input.lengths} {params.maxdistinter} {params.mindistinter} {params.maxgen} {params.mincirclen} {params.minleninter} {params.maxleninter} {output.igs}
        """


rule find_upstreams_coordinates_meta:
    input:
        rnaps = os.path.join(config['annotation_dir'], 'concatenated_rnaps_only.gff'),
        lengths = os.path.join(config['intergenic_dir'], 'chromosome_lengths.bed'),
        intergenics = os.path.join(config['intergenic_dir'], 'best_intergenics_meta.tsv')
    output:
        os.path.join(config['target_dir'], 'target_meta.bed')
    params:
        script = os.path.join(config['scripts'], 'get_upstream_bed_meta.py'),
        distance = config['max_dist_from_pol'],
        mincirclen = config['min_contig_len_circ']
    shell:
        """
        python {params.script} -r {input.rnaps} -i {input.intergenics} \
        -l {input.lengths} -b {params.distance} -o {output} -s {params.mincirclen}
        """


rule get_target_genes:
    input:
        bed = os.path.join(config['target_dir'], 'target.bed') if config['genomes_type'] == 'refseq'
        else os.path.join(config['target_dir'], 'target_meta.bed'),
        gff = os.path.join(config['annotation_dir'], 'concatenated.gff')
    output:
        gff = os.path.join(config['target_dir'], 'target.gff')
    conda: os.path.join(config['envs'], 'bedtools.yml')
    shell:
        """
        cat {input.bed} | cut -f 1-6 > temp_file.bed
        bedtools intersect -a {input.gff} -b temp_file.bed -s > {output}
        rm temp_file.bed
        """

rule get_target_faa:
    input:
        faa = os.path.join(config['annotation_dir'], 'concatenated.faa'),
        gff = os.path.join(config['target_dir'], 'target.gff'),
    output:
        faa_total = os.path.join(config['target_dir'], 'target.faa')
    params:
        script = os.path.join(config['scripts'], 'get_upstream_proteins_faa.py')
    conda: os.path.join(config['envs'], 'biopython.yml')
    shell:
        """
        python {params.script} -i {input.faa} -g {input.gff} -o {output.faa_total}
        """
