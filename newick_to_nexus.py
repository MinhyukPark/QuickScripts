import sys

import click
import dendropy


@click.command()
@click.option("--newick-path", required=True, type=click.Path(exists=True), help="Input newick file to be converted to Nexus format")
@click.option("--nexus-path", required=True, type=click.Path(exists=True), help="Output nexus file name")
def newick_to_nexus(newick_path, nexus_path):
    '''This program takes in a newick tree and converts it to a nexus format.
    '''
    tree = dendropy.Tree.get(path=newick_path, schema="newick")
    tree.write(path=nexus_path, schema="nexus")

if __name__ == "__main__":
    newick_to_nexus()
