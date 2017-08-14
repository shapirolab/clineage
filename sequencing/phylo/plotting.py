import copy
import contextlib
from sequencing.calling.queries.mutation_maps import filter_bipartition_loci, filter_mutation_map
from sequencing.calling.queries.formatters import get_cells_data_dict
from sequencing.phylo.csv_writers import write_cell_data_dict_to_file, print_mutation_dict_to_file
from functools import partial
from sequencing.phylo.matlab_wrappers import basic_tree_enrichment_and_plotting
from misc.utils import unlink, get_unique_path


@contextlib.contextmanager
def prep_plot_mutation_map_and_cell_data(srs, full_td, pl1, pl2, group_of_cell=lambda cell: cell.name):
    td = copy.deepcopy(full_td)
    td = filter_mutation_map(td, pl1, pl2)
    td = filter_bipartition_loci(td)
    with unlink(get_unique_path('tab')) as mutation_table_path:
        print_mutation_dict_to_file(td, mutation_table_path)
        plate_data_dict = get_cells_data_dict(srs, group_of_cell=group_of_cell)
        with unlink(get_unique_path('tab')) as cell_data_path:
            write_cell_data_dict_to_file(plate_data_dict, cell_data_path)
            plot = partial(
                basic_tree_enrichment_and_plotting,
                cell_data_path=cell_data_path,
                mutation_table_path=mutation_table_path,
                )
            yield mutation_table_path, plot
