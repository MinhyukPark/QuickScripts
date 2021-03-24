import sys

from Bio import SeqIO
import click

@click.command()
@click.option("--input-file", required=True, type=click.Path(exists=True), help="Input fasta file to be converted to phylip format")
@click.option("--output-file", required=True, type=click.Path(), help="Output phylip file name")
def fasta_to_phylip(input_file, output_file):
    '''This program takes in a fasta file and turns it into phylip file
    '''
    SeqIO.convert(input_file, "fasta", output_file, "phylip")

if __name__ == "__main__":
    fasta_to_phylip()
