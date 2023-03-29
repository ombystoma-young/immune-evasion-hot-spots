import os
import json

source_data = 'psi_blast_results/psi_blast_results.csv'
assemblies_dir = 'ncbi_dataset/data'
meta_dir = 'metadata'
stats_dir = 'stats'
pics_dir = 'pics'

genomes = os.listdir(assemblies_dir)
genomes.remove('assembly_data_report.jsonl')
genomes.remove('dataset_catalog.json')

os.makedirs(stats_dir, exist_ok=True)

rule all:
    input:
        stat=os.path.join(meta_dir,'stats_for_preso'),
        plot=os.path.join(pics_dir,'hists_gc_length.png'),
        drop=os.path.join(meta_dir,'genomes_to_drop.tsv')

# stats for each genome
rule obtain_stats:
    input:
        os.path.join(assemblies_dir, '{genome}', 'sequence_report.jsonl'),
        os.path.join(meta_dir, 'assembly_ids.txt')
    output:
        temp(os.path.join(stats_dir, '{genome}.stat'))
    threads: 8
    run:
        id_strain_dict = {}
        path_file = str(input[0])
        path_file_strain = str(input[1])
        path_out = str(output)

        with open(path_file_strain, 'r') as tsv_file:
            for line in tsv_file:
                comp = line.strip().split('\t')
                if len(comp) > 1:
                    id_strain_dict[comp[1]] = comp[0]

        with open(path_file, 'r') as json_file:
            js = json_file.readline()
            js = json.loads(js, object_hook=dict)

        aa_check = 'assemblyAccession' in js.keys()
        gcp_check = 'gcPercent' in js.keys()
        l_check = 'length' in js.keys()
        strain = id_strain_dict[js['assemblyAccession']]
        row = []
        if aa_check and gcp_check and l_check:
            row = [strain, str(js['assemblyAccession']), str(js['gcPercent']), str(js['length'])]
        elif not gcp_check and l_check:
            row = [strain, str(js['assemblyAccession']), 'None', str(js['length'])]
        elif not l_check and gcp_check:
            row = [strain, str(js['assemblyAccession']),  str(js['gcPercent']), 'None']
        else:
            row = [strain, str(js['assemblyAccession']), 'None', 'None']
        with open(path_out, 'w') as out_file:
            out_file.write('\t'.join(row))
            out_file.write('\n')

# collect stats for all genomes
rule collect_gc_length_stats:
    input:
        expand(os.path.join(stats_dir, '{genome}.stat'), genome=genomes)
    output:
        os.path.join(stats_dir, 'genomes_gc_length.statistics')
    threads: 8
    shell:
        """
        cat {input} > {output}
        """

# summary
rule dataset_stats_bash:
    input:
        stat=os.path.join(stats_dir, 'genomes_gc_length.statistics'),
        blast=source_data,
        meta=os.path.join(meta_dir, 'assembly_ids.txt')
    output:
        stat=os.path.join(meta_dir, 'stats_for_preso')
    shell:
        """
        sh scripts/dataset.stats.sh {input.stat} {input.blast} {input.meta} {output.stat}
        """
# hists
rule dataset_stats_plot:
    input:
        stat=os.path.join(stats_dir, 'genomes_gc_length.statistics')
    output:
        plot=os.path.join(pics_dir, 'hists_gc_length.png'),
        drop=os.path.join(meta_dir, 'genomes_to_drop.tsv')
    shell:
        """
        Rscript scripts/plot_initial_stats.R {input.stat} {output.plot} {output.drop}
        """