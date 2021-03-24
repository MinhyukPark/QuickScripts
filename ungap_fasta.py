from Bio import SeqIO
from Bio.Seq import Seq
import click


@click.command()
@click.option("--input-fasta", required=True, type=click.Path(exists=True), help="Input fasta file with gaps")
@click.option("--output-fasta", required=True, type=click.Path(writable=True), help="Output fasta file without gaps")
def ungap_fasta(input_fasta, output_fasta):
    '''This program takes in a fasta file and removes all gaps
    '''
    with open(output_fasta, "w") as f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            record.seq = record.seq.ungap("-")
            SeqIO.write(record, f, "fasta")

if __name__ == "__main__":
    ungap_fasta()
