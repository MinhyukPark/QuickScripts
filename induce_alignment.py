import random
import sys

from Bio import SeqIO
import click

@click.command()
@click.option("--input-sequences", required=True, type=click.Path(exists=True), help="Input fasta sequence file")
@click.option("--to-keep-sequences", required=True, type=click.Path(exists=True), help="Input fasta sequence file to keep")
@click.option("--output-prefix", required=True, type=click.Path(), help="Ouput sequence file prefix")
def induce_alignment(input_sequences, to_keep_sequences, output_prefix):
    """This script induces the input alignment on a sequence file
    """
    to_keep_set = {}
    for to_keep_sequence in SeqIO.parse(to_keep_sequences, "fasta"):
        to_keep_set[to_keep_sequence.id] = 1
    with open(output_prefix + ".keep", "w") as keep_f:
        with open(output_prefix + ".discard", "w") as discard_f:
            for sequence in SeqIO.parse(input_sequences, "fasta"):
                if(sequence.id in to_keep_set):
                    SeqIO.write(sequence, keep_f, "fasta")
                else:
                    SeqIO.write(sequence, discard_f, "fasta")

if __name__ == "__main__":
    induce_alignment()
