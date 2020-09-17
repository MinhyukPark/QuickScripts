import click
import dendropy

@click.command()
@click.option('--input-tree', required=True, type=click.Path(exists=True))
@click.option('--output-tree', required=True)
def remove_labels(input_tree, output_tree):
    tree = dendropy.Tree.get(path=input_tree, schema="newick")
    tree.write(path=output_tree, schema="newick", suppress_internal_taxon_labels=True, suppress_internal_node_labels=True, suppress_edge_lengths=True)

if __name__ == '__main__':
    remove_labels()
