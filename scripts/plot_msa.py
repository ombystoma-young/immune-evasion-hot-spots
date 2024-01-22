import re
import argparse

from Bio.SeqIO.FastaIO import SimpleFastaParser
from pymsaviz import MsaViz


def find_gaps(msa_file: str, seq_id: str, size=5) -> list:
    """
    finds coordinates of gaps
    :param msa_file: (str) path to msa file
    :param seq_id: (str) seq_id, where is needed to find gaps
    :param size: minimal size of a gap, 10
    :return: (list) of (tuples), corresponding to start and end of gap regions
    """
    gaps = []
    with open(msa_file, 'rt') as fasta:
        for record in SimpleFastaParser(fasta):
            if record[0] == seq_id:
                pattern = '-{' + str(size) + ',}'
                for gap in re.finditer(pattern, record[1]):
                    gap_coord_1based = gap.span()[0] + 1, gap.span()[1]
                    gaps.append(gap_coord_1based)
    print(gaps)
    return gaps


def plot_msa(msa_file: str, out_file: str, insertions: list = None, regions: list = None,
             mark_positions: list = None) -> None:
    """

    :param msa_file: (str) path to file with MSA
    :param insertions: (list) of (tuples), corresponding to start and end of insertion regions
    :param regions: (dict) of (tuples), corresponding to start and end of insertion regions
    :param out_file: (str) path to output figure
    :return: None
    """
    mv = MsaViz(msa_file, wrap_length=46, show_count=True, color_scheme='Identity',
                label_type='description', )
    mv.set_plot_params(identity_color='#7570B3')
    if mark_positions is not None:
        mv.add_markers(positions=mark_positions, color='#d95f02')
    if insertions is not None:
        for insertion in insertions:
            mv.add_text_annotation(insertion, "Insert Region", text_color="#E7298A", range_color="#E7298A")
    if regions is not None:
        for region_name in regions.keys():
            mv.add_text_annotation(regions[region_name][0], region_name, text_color=regions[region_name][1], range_color=regions[region_name][1])
    mv.savefig(out_file, dpi=300)


def parse_args():
    parser = argparse.ArgumentParser('Plot alignment')
    parser.add_argument('-i', '--input', default=None, nargs='?',
                        help='input msa file')
    parser.add_argument('-o', '--output', default=None, nargs='?')
    return parser.parse_args()


if __name__ == '__main__':
    msa_file = parse_args().input
    output_path = parse_args().output
    if 'samase' in msa_file:
        markers = [73, 76, 87, 113]
    else:
        markers = None
    plot_msa(msa_file, insertions=None, out_file=output_path, regions=None, mark_positions=markers)