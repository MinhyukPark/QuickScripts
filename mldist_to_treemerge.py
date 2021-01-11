import click

@click.command()
@click.option("--input-filename", required=True, type=click.Path(exists=True), help="Input iqtree mldist file")
@click.option("--output-path", required=True, type=click.Path(), help="Output path for taxa list")
def mldist_to_treemerge(input_filename, output_path):
    with open(input_filename, "r") as f:
        with open(output_path, "w") as fw:
            input_line_counter = 0
            for line in f:
                if(input_line_counter > 0):
                    fw.write(line.split()[0])
                    fw.write("\n")
                input_line_counter += 1

if __name__ == "__main__":
    mldist_to_treemerge();
