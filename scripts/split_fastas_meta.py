import os
import argparse

from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_meta_annotate_genomes" pipeline. '
                                                 'Split multifasta file into a small fasta files '
                                                 'where one file corresponds to one genome. '
                                                 'If specified, performs the filtering of contig entries,'
                                                 ' based on list file with contigs passed the filter.')
    parser.add_argument('-i', '--input', default=None, type=str, nargs='?',
                        help='path to input fasta file')
    parser.add_argument('-f', '--filter', default=None, type=str, nargs='?',
                        help='path to file with desired contigs. (default: None, no filtering provided)')
    parser.add_argument('-b', '--batches', default=1, type=int, nargs='?',
                        help='Number of batches (output files)')
    parser.add_argument('-l', '--minlen', default=1, type=int, nargs='?',
                        help='Minimal length of contig')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?',
                        help='path to output directory, which will contains files with ".fna" extension'
                             ' and filenames corresponding to output contig names')
    return parser.parse_args()


# def batch_iterator_no_filtering(iterator, batch_size: int, minlen: int):
#     """Returns lists of length batch_size.
#
#     This can be used on any iterator, for example to batch up
#     SeqRecord objects from Bio.SeqIO.parse(...), or to batch
#     Alignment objects from Bio.Align.parse(...), or simply
#     lines from a file handle.
#
#     This is a generator function, and it returns lists of the
#     entries from the supplied iterator.  Each list will have
#     batch_size entries, although the final list may be shorter.
#     SOURCE: https://biopython.org/wiki/Split_large_file
#     """
#     batch = []
#     for entry in iterator:
#         print(entry)
#         if len(entry) > minlen:
#             batch.append(entry)
#             if len(batch) == batch_size:
#                 yield batch
#                 batch = []
#     if batch:
#         yield batch


def batch_iterator(iterator, batch_size: int, minlen: int, filtering: list | None):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.Align.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    SOURCE: https://biopython.org/wiki/Split_large_file
    """
    batch = []
    if filtering is not None:
        for entry in iterator:
            if len(entry.seq) > minlen:
                if entry.id in filtering:
                    batch.append(entry)
                    if len(batch) == batch_size:
                        yield batch
                        batch = []
        if batch:
            yield batch
    else:
        for entry in iterator:
            if len(entry.seq) > minlen:
                batch.append(entry)
                if len(batch) == batch_size:
                    yield batch
                    batch = []


def split_fasta(in_path: str, num_batches: int, min_len: int, out_dir: str, filtering: list | None) -> None:
    """
    biopython SeqIO based function which parse input fasta file and writes FILTERED contigs into a
    separated files in desired directory
    :param num_batches: (int) number of batches
    :param in_path: (str) path to input file in fasta format
    :param out_dir: (str) path to output directory which will contain output files
    :param filtering: (list) of contigs which should be written into an output directory
    (in a particular pipeline obtainted from CheckV completeness check)
    :return: None
    """
    os.makedirs(out_dir, exist_ok=True)
    batch_size = len(list(SeqIO.parse(in_path, "fasta"))) // (num_batches - 1)
    record_iter = SeqIO.parse(in_path, "fasta")
    for i, batch in enumerate(batch_iterator(record_iter, batch_size, min_len, filtering)):
        out_path = os.path.join(out_dir, f"batch_{i + 1}.fna")
        with open(out_path, "wt") as output_handle:
            SeqIO.write(batch, output_handle, "fasta")


def read_txt_file(in_path: str | None) -> list | None:
    """
    reads txt file into a list of entries
    :param in_path: (str) path of input file
    :return: (list) of entries from file
    """
    if in_path is None:
        return None
    entries = []
    with open(in_path, 'rt') as in_f:
        for line in in_f:
            entries.append(line.strip())
    return entries


if __name__ == '__main__':
    in_path = parse_args().input
    batches_num = parse_args().batches
    min_len = parse_args().minlen
    out_path = parse_args().output
    passed_contigs_path = parse_args().filter

    passed_contigs = read_txt_file(in_path=passed_contigs_path)
    split_fasta(in_path=in_path,
                min_len=min_len,
                num_batches=batches_num,
                out_dir=out_path,
                filtering=passed_contigs)
