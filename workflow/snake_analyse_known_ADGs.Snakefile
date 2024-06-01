import os

configfile: 'config_autographiviridae.yaml'

# CREATE FOLDERS
os.makedirs(config['known_interest_dir'], exist_ok=True)
os.makedirs(config['alignments_dir'], exist_ok=True)
os.makedirs(config['alignments_trimmed_dir'], exist_ok=True)
os.makedirs(config['trees_dir'], exist_ok=True)


def create_search_string(nums: list) -> str:
    nums_list = [f'clu={num};' for num in nums]
    nums_str = '|'.join(nums_list)
    return nums_str



rnaps = [1264, 15258]
samases = [3643, 9717, 12011, 17333, 17415, 7377,
           175, 8990, 12408, 4202, 10755, 2290,
           15688, 16209, 1511]
ocrs = [14462, 15652, 17328, 6718,
        6958, 8013, 8030, 8611]
kinases = [
    9668, 10265, 3025, 10700,
    2969,
    12312, 15065, 1881, 13989
          ]

mreb_inhibitors = [
                    1293, 1726, 17321, 4442,
                    54, 6799, 7355, 7429, 7503, 838
                    ]



clusters = ['ocr', 'samase', 'rnap', 'arda']
clu_nums_ocr_str = create_search_string(ocrs)
clu_nums_rnaps_str = create_search_string(rnaps)
clu_nums_samase_str = create_search_string(samases)

clu_nums_kinases_str = create_search_string(kinases)
clu_nums_06_str = create_search_string(mreb_inhibitors)

clu_num = {'ocr': clu_nums_ocr_str,
           'samase': clu_nums_samase_str,
           'rnap': clu_nums_rnaps_str,
           '06': clu_nums_06_str,
           'kinase': clu_nums_kinases_str
           }

clu_files = {'ocr': ocrs,
           'samase': samases,
           'rnap': rnaps,
           '06': mreb_inhibitors,
           'kinase': kinases
           #  'interest': [3030, 3020, 1677, 516, 3704]
             }

clusters_exp = ['ocr', 'samase', 'rnap', 'kinase', '06']
rule all:
    input:
        expand([os.path.join(config['trees_dir'], '{family}_fasttree.treefile'),
                os.path.join(config['alignments_trimmed_dir'], 'trimmed_{family}.mafft.faa')],
                family=['rnap']),
        expand(os.path.join(config['alignments_dir'], 'clu_reprs_mafft_{family}.svg'),
             family=clusters_exp)

rule select_cluster_representatives:
    input:
        os.path.join(config['target_dir'], 'target_with_clusters_phrogs.gff')
    output:
        os.path.join(config['known_interest_dir'], 'target_{family}.gff')
    params:
        keyword=lambda wildcards: clu_num[f'{wildcards.family}']
    shell:
        """
        cat {input} | grep -P "{params.keyword}" > {output} 
        """


rule get_protein_ids:
    input:
        os.path.join(config['similarity_dir'], 'target_with_clusters_phrogs_within_communities.gff')
    output:
        os.path.join(config['known_interest_dir'], 'protein_ids_{family}.tsv')
    shell:
        """
        cat {input} | tr -s "\t" ";" | cut -f 1,9 -d ";" | tr -s "=" ";" | cut -f 1,3 -d ";"  | tr -s ";" "\t" > {output}
        """


# BLOCK ALIGN ANTIDEFENCE PROTEINS FOR BEST DATASET: MAFFT + TRIMAL + IQTREE
rule get_cluster_faa:
    input:
        faa = os.path.join(config['target_dir'], 'target.faa'),
        tsv = os.path.join(config['known_interest_dir'], 'protein_ids_{family}.tsv')
    output:
        faa = os.path.join(config['known_interest_dir'],  'upsteam_{family}.faa')
    params:
        script = os.path.join(config['scripts'], 'get_cluster_faa.py')
    shell:
        """
       python {params.script} {input.faa} {input.tsv} {output} 
       """


rule align_proteins:
    input:
        os.path.join(config['known_interest_dir'],  'upsteam_{family}.faa')
    output:
        os.path.join(config['alignments_dir'],  'mafft_{family}.faa')
    conda: os.path.join(config['envs'], 'mafft.yml')
    threads: config['maxthreads']
    shell:
        """
         mafft --thread {threads} --maxiterate 1000 --globalpair {input} > {output} 
        """


rule filter_alignment:
    input:
        os.path.join(config['alignments_dir'],  'mafft_{family}.faa')
    output:
        os.path.join(config['alignments_trimmed_dir'], 'trimmed_{family}.mafft.faa')
    conda: os.path.join(config['envs'], 'trimal.yml')
    shell:
        """
        trimal -in {input} -out {output} -automated1
        """


rule build_tree_fast_tree:
    input:
        os.path.join(config['alignments_trimmed_dir'], 'trimmed_{family}.mafft.faa')
    output:
        os.path.join(config['trees_dir'], '{family}_fasttree.treefile')
    conda:
        os.path.join(config['envs'], 'fasttree.yml')
    shell:
        """
        FastTree {input} > {output}
        """

rule build_tree_iqtree:
    input:
        os.path.join(config['alignments_trimmed_dir'], 'trimmed_{family}.mafft.faa')
    output:
        os.path.join(config['trees_dir'], '{family}_bootstrap_model_selection.iqtree.treefile')
    params:
        pref=os.path.join(config['trees_dir'], '{family}_bootstrap_model_selection.iqtree'),
        model='MFP',
        ubootstrap=10000
    threads: config['maxthreads']
    conda: os.path.join(config['envs'], 'iqtree2.yml')
    shell:
        """
        iqtree2 -T AUTO -m {params.model} -s {input[0]} --prefix {params.pref} -B {params.ubootstrap} -redo
        """

rule extract_clu_representatives:
    input:
        fa = config['target_msa_unpacked_dir'],
        source_fa = os.path.join(config['target_dir'], 'target.faa')
    output:
        os.path.join(config['known_interest_dir'], 'clu_reprs_{family}.faa')
    params: clu_nums = lambda wc: ','.join([str(i) for i in clu_files[wc.family]]),
            script = os.path.join(config['scripts'], 'extract_clu_representatives.py')
    shell:
        """
        python3 {params.script} -d {input.fa} -f {input.source_fa} -s {params.clu_nums} -o {output}
        """
#
rule align_representatives:
    input:
        os.path.join(config['known_interest_dir'], 'clu_reprs_{family}.faa')
    output:
        os.path.join(config['alignments_dir'], 'clu_reprs_mafft_{family}.faa')
    conda: os.path.join(config['envs'], 'mafft.yml')
    threads: 10
    shell:
        """
         mafft --thread {threads} --maxiterate 1000 --localpair {input} > {output}
        """

rule plot_alignment:
    input:
        os.path.join(config['alignments_dir'], 'clu_reprs_mafft_{family}.faa')
    output:
        os.path.join(config['alignments_dir'], 'clu_reprs_mafft_{family}.svg')
    params:
        script = os.path.join(config['scripts'], 'plot_msa.py')
    shell:
        """
        python {params.script} -i {input} -o {output}
        """