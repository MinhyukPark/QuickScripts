import random
import sys

from Bio import SeqIO
import click

@click.command()
@click.option("--input-sequence", required=True, type=click.Path(exists=True), help="Input fasta sequence file")
@click.option("--output-sequence", required=True, type=click.Path(), help="Ouput sequence file name")
@click.option("--num-sequence", required=True, type=int, help="Number of sequences to randomly sample")
def random_sample_sequence(input_sequence, output_sequence, num_sequence):
    """This script induces the input tree on the sequence file
    """
    sequences = list(SeqIO.parse(open(input_sequence), "fasta"))
    random.shuffle(sequences)
    SeqIO.write(sequences[:num_sequence], open(output_sequence, "w"), "fasta")


if __name__ == "__main__":
    random_sample_sequence()
