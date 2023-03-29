import os


def read_gff(path_to_gff):
    cds_ids = []
    with open(path_to_gff, 'r') as gff:
        for line in gff:
            cds = line.strip().split('\t')
            attribute = cds[-1]
            attribute_list = attribute.split(';')
            cds_id = attribute_list[0].split('=')[1]
            cds_ids.append(cds_id)
    return set(cds_ids)


def read_fa(gff_p, faa_path):
    cds_ids = read_gff(gff_p)
    sequences = {}
    name = None
    with open(faa_path, 'r') as in_faa:
        for line in in_faa:
            if line.startswith('>'):
                if name in cds_ids:
                    sequences[name] = seq
                name = line.strip().split(' ')[0][1:]
                seq = ''
            else:
                seq += line.strip()
        if name in cds_ids:
            sequences[name] = seq
    return sequences


def find_upstream_faa(gff_path, faa_input_path, faa_output_path):
    seqs = read_fa(gff_path, faa_input_path)
    with open(faa_output_path, 'w') as out_faa:
        for cds_id in seqs.keys():
            out_faa.write(f'>{cds_id}\n')
            out_faa.write(f'{seqs[cds_id]}\n')


if __name__ == '__main__':
    folder = 'search_upstream'
    find_upstream_faa(gff_path=os.path.join(folder, 'upstream.gff'),
                      faa_input_path=os.path.join(folder,'all_genomes.faa'),
                      faa_output_path=os.path.join(folder, 'upstream.faa'))