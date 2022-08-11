import click
import dendropy

@click.command()
@click.option("--input-filename", required=True, type=click.Path(exists=True), help="Input newick tree file")
@click.option("--output-path", required=True, type=click.Path(), help="Output path for distance matrix")
@click.option("--dist-type", type=click.Choice(["node", "brlen"], case_sensitive=False), required=True, help="Type of distance to calculate")
def tree_to_dist(input_filename, output_path, dist_type):
    tree_to_dist_helper(input_filename, output_path, dist_type)

def tree_to_dist_helper(input_filename, output_path, dist_type):
    tree = dendropy.Tree.get(path=input_filename, schema="newick")

    is_weighted_edge_distances=None
    if(dist_type == "node"):
        is_weighted_edge_distances=False
    elif(dist_type == "brlen"):
        is_weighted_edge_distances=True

    pdc = tree.phylogenetic_distance_matrix()
    with open(output_path, "w") as f:
        f.write(str(len(tree.leaf_nodes())))
        f.write("\n")
        for row_taxon in tree.taxon_namespace:
            f.write(row_taxon.label + " ")
            for col_taxon in tree.taxon_namespace:
                f.write(str(pdc.distance(row_taxon, col_taxon, is_weighted_edge_distances)) + " ")
            f.write("\n")

if __name__ == "__main__":
    tree_to_dist();
