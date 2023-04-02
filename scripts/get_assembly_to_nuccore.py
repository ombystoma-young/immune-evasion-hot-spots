import os


def read_jsonl(path: str) -> list:
    entries = []
    with open(path, 'r') as jsonl:
        for line in jsonl:
            entry = line.strip()[1:-1].split(',')
            entry_dict = {el.split(':')[0]: el.split(':')[1] for el in entry}
            entries.append(entry_dict)
    return entries


def get_chr_length_couples(entries: list) -> list:
    chr_assemblies = []
    for entry in entries:
        assembly = entry['"assemblyAccession"']
        if '"refseqAccession"' in entry.keys():
            chrm_refseq = entry['"refseqAccession"']
            chr_assembly = (chrm_refseq, assembly)
            # chr_assemblies.append(chr_assembly)
        else:
            chrm_genbank = entry['"genbankAccession"']
            chr_assembly = (chrm_genbank, assembly)
        chr_assemblies.append(chr_assembly)
    return chr_assemblies


def write_tsv(couples_list: list, tsv_name: str) -> None:
    with open(tsv_name, 'w') as bed_file:
        for couple in couples_list:
            nuccore = couple[0][1:-1]
            assembly = couple[1][1:-1]
            line = "\t".join([nuccore, assembly])
            bed_file.write(line)
            bed_file.write('\n')


if __name__ == '__main__':
    in_jsonl = os.path.join('promoters_search', 'all_sequence_report.jsonl')
    out_tsv = os.path.join('metadata', 'assembly_nuccore.tsv')
    e = read_jsonl(in_jsonl)
    couples = get_chr_length_couples(e)
    write_tsv(couples_list=couples, tsv_name=out_tsv)
