import sys
sys.path.append('/home/dcsoft/s/Ofir/igor_tree_reconstruction_20170731/')
from importalg import importalg
from backgraund import backgraund7
from read_data import parse_mutations_table, get_cells_and_root, parse_distance_matrix, write_distance_matrix


def calculate_igors_tree(mutation_table_path, output_path, dist="ig1", eps=0.001):
    calling = parse_mutations_table(mutation_table_path, inverse=True)
    root, all_cells = get_cells_and_root(calling)
    cells = list(all_cells)
    matrix = backgraund7(cells, dist)
    newick1, list2g = importalg(matrix, eps)
    if root is not None:
        newick = "({})'{}':0;".format(newick1, root)
    else:
        newick = "{};".format(newick1)
    with open(output_path, 'w') as f:
        f.write(newick)
