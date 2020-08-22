import sys

import dendropy

input_path = sys.argv[1]
output_path = sys.argv[2]

tree = dendropy.Tree.get(path=input_path, schema="newick")
tree.write(path=output_path, schema="nexus")
