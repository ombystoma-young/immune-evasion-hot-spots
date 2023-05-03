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
    chr_lengths = []
    for entry in entries:
        if '"refseqAccession"' in entry.keys():
            chrm = entry['"refseqAccession"']
        else:
            chrm = entry['"genbankAccession"']
        length = entry['"length"']
        chr_len = (chrm, length)
        chr_lengths.append(chr_len)
    return chr_lengths


def write_bed(couples_list: list, bed_name: str) -> None:
    with open(bed_name, 'w') as bed_file:
        for couple in couples_list:
            chrm = couple[0][1:-1]
            length = couple[1]
            line = "\t".join([chrm, str(0), length])
            bed_file.write(line)
            bed_file.write('\n')


if __name__ == '__main__':
    in_jsonl = os.path.join('promoters_search', 'all_sequence_report.jsonl')
    out_bed = os.path.join('promoters_search', 'chromosome_lengths.bed')
    e = read_jsonl(in_jsonl)
    couples = get_chr_length_couples(e)
    write_bed(couples_list=couples, bed_name=out_bed)
