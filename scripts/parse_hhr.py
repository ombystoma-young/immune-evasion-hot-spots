import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_process_clustering_results" pipeline. '
                                                 'Parse hhr file into a table')
    parser.add_argument('-i', '--input', default=None, type=str, nargs='?',
                        help='path to input hhr file')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output table file')
    return parser.parse_args()


def parse_hhr_file(file_path):
    hits = []
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                # Parse query name
                query = line.strip().split()[1]
                query_name = query.split('|')[0].split('-')[1]
            elif i == 8:
                # Parse header information
                header = line.strip().split()[:-4] + ['Query HMM', 'Template HMM', 'Length']
            elif i > 7:
                if len(line.strip()) == 0:
                    break
                else:
                    # Parse hit information
                    entry = line.strip().split()
                    if entry[2].startswith('n='):
                        entry.pop(2)
                    if entry[2].startswith('D'):
                        entry.pop(2)
                    if len(entry) < 11:
                        if '(' in entry[-1]:
                            problematic_range = entry[-1]
                            entry[-1] = problematic_range.split('(')[0]
                            entry.append('(' + problematic_range.split('(')[1])
                        else:
                            raise ValueError(f'{query_name, line}: wrong format')
                    elif len(entry) > 11:
                        raise ValueError(f'{query_name, entry}: wrong format')

                    entry[1] = entry[1].split('|')[0].split('-')[1]
                    hits.append(entry)

    # Create pandas dataframe
    df = pd.DataFrame(hits, columns=header)
    df['Query HMM'] = df['Query HMM'].apply(lambda x: x.split('-'))
    df['Query_start'] = df['Query HMM'].apply(lambda x: int(x[0]))
    df['Query_end'] = df['Query HMM'].apply(lambda x: int(x[1]))
    df.drop(columns=['Query HMM'], inplace=True)
    df['Template HMM'] = df['Template HMM'].apply(lambda x: x.split('(')[0])
    df['Template_start'] = df['Template HMM'].apply(lambda x: int(x.split('-')[0]))
    df['Template_end'] = df['Template HMM'].apply(lambda x: int(x.split('-')[1]))
    df.drop(columns=['Template HMM'], inplace=True)
    # Add header information as metadata
    df['Query'] = query_name
    return df


if __name__ == '__main__':
    in_file = parse_args().input
    out_file = parse_args().output
    # in_file = '../data/clans_pre_info/635.hhr'
    hh_results = parse_hhr_file(in_file)
    # print(hh_results)
    hh_results.to_csv(out_file, sep='\t', index=False)