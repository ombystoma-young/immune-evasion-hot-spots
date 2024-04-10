import os

configfile: 'config_autographiviridae_refseq.yaml'

# CREATE FOLDERS
os.makedirs(config['known_interest_dir'], exist_ok=True)
os.makedirs(config['alignments_dir'], exist_ok=True)
os.makedirs(config['alignments_trimmed_dir'], exist_ok=True)
os.makedirs(config['trees_dir'], exist_ok=True)


def create_search_string(nums: list) -> str:
    nums_list = [f'clu={num};' for num in nums]
    nums_str = '|'.join(nums_list)
    return nums_str


#rnaps = [3832, 3398, 2681, 2172, 803, 156]
#samases = [392, 516, 2538, 3167, 3357, 3486, 3509, 3704, 1677]
#ocrs = [3814, 3136, 3132, 3030, 3020, 1911, 1677]
#kinases = [21, 75, 237, 416, 622, 1390, 1458,
#            1659, 2684, 2865, 3132, 3521, 3566, 3654]

rnaps = [7905, 3093, 4851, 6262, 2689, 1423, 4636, 5335, 3763]
samases = [4422, 7327, 3809]
ocrs = [3892, 4935, 7170, 1528]
kinases = [4877, 364, 3817, 4632, 2103, 3212]
arda = [8435, 3512, 6059]

#rnaps = [665, 2336, 499, 2679, 2003, 3009, 285, 558, 1244, 1175, 2670, 382, 677]
#samases = [365, 2059, 957]
#ocrs = [851, 19]
#arda = [357, 2015, 1747]


#clusters = ['ocr', 'samase', 'rnap', 'arda', 'kinase']
clusters = ['ocr', 'samase', 'rnap', 'arda']
clu_nums_ocr_str = create_search_string(ocrs)
clu_nums_rnaps_str = create_search_string(rnaps)
clu_nums_samase_str = create_search_string(samases)
clu_nums_arda_str = create_search_string(arda)
clu_nums_kinases_str = create_search_string(kinases)

clu_num = {'ocr': clu_nums_ocr_str,
           'samase': clu_nums_samase_str,
           'rnap': clu_nums_rnaps_str,
           'arda': clu_nums_arda_str,
           'kinase': clu_nums_kinases_str
           }

clu_files = {'ocr': ocrs,
           'samase': samases,
           'rnap': rnaps,
           'arda': arda,
           'kinase': kinases
           #  'interest': [3030, 3020, 1677, 516, 3704]
             }
# clusters_exp = ['ocr', 'samase', 'rnap', 'kinase', 'interest']
clusters_exp = ['ocr', 'samase', 'rnap', 'arda']
rule all:
    input:
        expand([os.path.join(config['trees_dir'], '{family}_fasttree.treefile'),
                os.path.join(config['alignments_trimmed_dir'], 'trimmed_{family}.mafft.faa')],
                family=clusters),
        # expand(os.path.join(config['alignments_dir'],'representatives_mafft_{family}.svg'),
        #     family=clusters_exp)

rule select_cluster_representatives:
    input:
        os.path.join(config['upstreams_dir'], 'early_with_clusters_phrogs.gff')
    output:
        os.path.join(config['known_interest_dir'], 'upstreams_{family}.gff')
    params:
        keyword=lambda wildcards: clu_num[f'{wildcards.family}']
    shell:
        """
        cat {input} | grep -P "{params.keyword}" > {output} 
        """


rule get_protein_ids:
    input:
        os.path.join(config['similarity_dir'], 'early_with_clusters_phrogs_within_communities.gff')
    output:
        os.path.join(config['known_interest_dir'], 'protein_ids_{family}.tsv')
    shell:
        """
        cat {input} | tr -s "\t" ";" | cut -f 1,9 -d ";" | tr -s "=" ";" | cut -f 1,3 -d ";"  | tr -s ";" "\t" > {output}
        """


# BLOCK ALIGN ANTIDEFENCE PROTEINS FOR BEST DATASET: MAFFT + TRIMAL + IQTREE
rule get_cluster_faa:
    input:
        faa = os.path.join(config['upstreams_dir'], 'early.faa'),
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

rule extract_representatives:
    input:
        fa = config['early_msa_unpacked_dir']
    output:
        os.path.join(config['known_interest_dir'], 'representatives_{family}.faa')
    params: clu_nums = lambda wc: ','.join([str(i) for i in clu_files[wc.family]]),
            script = os.path.join(config['scripts'], 'extract_clu_representatives.py')
    shell:
        """
        python3 {params.script} -d {input} -s {params.clu_nums} -o {output}
        """
#
rule align_representatives:
    input:
        os.path.join(config['known_interest_dir'], 'representatives_{family}.faa')
    output:
        os.path.join(config['alignments_dir'], 'representatives_mafft_{family}.faa')
    conda: os.path.join(config['envs'], 'mafft.yml')
    threads: 10
    shell:
        """
         mafft --thread {threads} --maxiterate 1000 --localpair {input} > {output}
        """

rule plot_alignment:
    input:
        os.path.join(config['alignments_dir'], 'representatives_mafft_{family}.faa')
    output:
        os.path.join(config['alignments_dir'], 'representatives_mafft_{family}.svg')
    params:
        script = os.path.join(config['scripts'], 'plot_msa.py')
    shell:
        """
        python {params.script} -i {input} -o {output}
        """