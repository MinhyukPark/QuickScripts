import sys

import click


def relabel_distance_matrix_helper(input_matrix, tax_list, output_matrix):
    tax_map = {}
    with open(tax_list, "r") as f:
        line_counter = 0
        for line in f:
            tax_map[str(line_counter)] = line.strip()
            line_counter += 1
    print(tax_map)
    with open(input_matrix, "r") as f:
        with open(output_matrix, "w") as fw:
            line_counter = 0
            for line in f:
                if line_counter == 0:
                    fw.write(line + "\n")
                else:
                    current_line_arr = line.split()
                    current_line_arr[0] = tax_map[str(line_counter - 1)]
                    current_line = " ".join(current_line_arr)
                    fw.write(current_line + "\n")
                line_counter += 1

@click.command()
@click.option("--input-matrix", required=True, type=click.Path(exists=True), help="The input tree file in newick format")
@click.option("--tax-list", required=True, type=click.Path(exists=True), help="Taxlist for the input tree")
@click.option("--output-matrix", required=True, type=str, help="Output file path for the relabeled matrix")
def relabel_distance_matrix(input_matrix, tax_list, output_matrix):
    """This script relabel the matrix according to the input taxlist
    """
    relabel_distance_matrix_helper(input_matrix, tax_list, output_matrix)


if __name__ == "__main__":
    relabel_distance_matrix()
