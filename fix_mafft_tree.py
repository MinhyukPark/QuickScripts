import click
import dendropy

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="The input tree file in newick format")
@click.option("--mafft-output", required=True, type=click.Path(exists=True), help="The output file from running mafft")
@click.option("--output-file", required=True, type=str, help="Output file path for the induced subtree")
def join_multiline_tree(input_tree, mafft_output, output_file):
    """This script joins a multiline newick tree and adds a semicolon at the end because MAFFT and also relabels the trees
    """
    taxa_mapping = {}
    taxa_count = 1
    with open(mafft_output, "r") as mafft_in:
        for line in mafft_in:
            if line[0] == ">":
                taxa_label = line[1:].strip()
                taxa_mapping[taxa_count] = taxa_label
                taxa_count += 1

    with open(input_tree, "r") as f_in:
        with open(output_file + ".wrong_labels", "w") as f_out:
            for line in f_in:
                f_out.write(line.strip())
            f_out.write(";");

    mafft_tree = dendropy.Tree.get(path=output_file + ".wrong_labels", schema="newick")
    for node in mafft_tree.leaf_nodes():
        node.taxon.label = taxa_mapping[int(node.taxon.label)]
    mafft_tree.write(path=output_file, schema="newick")

if __name__ == "__main__":
    join_multiline_tree()
