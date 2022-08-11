import click

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="The input tree file in newick format")
@click.option("--output-file", required=True, type=str, help="Output file path for the induced subtree")
def join_multiline_tree(input_tree, output_file):
    """This script joins a multiline newick tree and adds a semicolon at the end because MAFFT
    """
    with open(input_tree, "r") as f_in:
        with open(output_file, "w") as f_out:
            for line in f_in:
                f_out.write(line.strip())
            f_out.write(";");

if __name__ == "__main__":
    join_multiline_tree()
