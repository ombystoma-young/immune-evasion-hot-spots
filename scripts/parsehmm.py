import pandas as pd
import argparse


def parse_hmm(hmm_data: str) -> list:
    lines = []
    with open(hmm_data, 'r') as table:
        for line in table:
            if not line.startswith('#'):
                lines.append(line.strip().split())
    entries = []     
    for line in lines:
        entry = line[:18]
        description = ' '.join(line[18:])
        entry.append(description)
        entries.append(entry)
    return entries


def convert_to_readable_pd(hmm_entries: list):
    col_names =  ['target_name', 'target_accession',
              'query_name',  'query_accession',
              'full_seq_E_value', 'full_seq_score', 'full_seq_bias',
              'best_one_domain_E_value', 'best_one_domain_score', 'best_one_domain_bias',
              'dom_num_est_exp', 'dom_num_est_reg',
              'dom_num_est_clu', 'dom_num_est_ov',
              'dom_num_est_env', 'dom_num_est_dom',
              'dom_num_est_rep', 'dom_num_est_inc',
              'description_of_target']
    cols_to_drop = ['target_accession',               
                'dom_num_est_exp',
                'dom_num_est_reg',
                'dom_num_est_clu',
                'dom_num_est_ov',
                'dom_num_est_env',
                'dom_num_est_dom',
                'dom_num_est_rep',
                'dom_num_est_inc', 
                'full_seq_score',
                'full_seq_bias',
                'best_one_domain_score',
                'best_one_domain_bias']
    
    data = pd.DataFrame(hmm_entries, columns=col_names)
    data = data.drop(columns=cols_to_drop)
    data.full_seq_E_value = pd.to_numeric(data.full_seq_E_value)
    data.best_one_domain_E_value = pd.to_numeric(data.best_one_domain_E_value)
    return data
    

def extract_best_hits(data):
    data_filtered = data.sort_values('full_seq_E_value', ascending=True).drop_duplicates(['target_name'])
    return data_filtered
    
    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', default=None, nargs='?')
    parser.add_argument('--extractbest', default=False)
    parser.add_argument('--outfile', default=None, nargs='?')
    return parser.parse_args()
    
    
if __name__ == '__main__':
    infile = parse_args().infile
    extractbest = parse_args().extractbest
    outfile = parse_args().outfile
    
    hmms = parse_hmm(hmm_data=infile)
    data = convert_to_readable_pd(hmms)
    if extractbest:
        data = extract_best_hits(data)
    data.to_csv(outfile, sep='\t', index=False)
