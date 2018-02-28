import copy
import contextlib
from sequencing.calling.queries.mutation_maps import filter_bipartition_loci, filter_mutation_map
from sequencing.calling.queries.formatters import get_cells_data_dict
from sequencing.phylo.csv_writers import write_cell_data_dict_to_file, print_mutation_dict_to_file
from functools import partial
from sequencing.phylo.matlab_wrappers import basic_tree_enrichment_and_plotting
from misc.utils import unlink, relaxed_unlink, get_unique_path


@contextlib.contextmanager
def prep_mutation_map_and_cell_data(srs, full_td, pl1, pl2, group_of_cell=lambda cell: cell.name):
    """
    The function
    :param srs: List of SampleReads
    :param full_td: Dict of cells with their loci length
    :param pl1:
    :param pl2:
    :param group_of_cell: list of the cell groups
    :return: (yield) td -
            mutation_table_path -
            cell_data_path -
            plot -
            output_plot -
    """
    td = copy.deepcopy(full_td)
    td = filter_mutation_map(td, pl1, pl2)
    td = filter_bipartition_loci(td)
    with relaxed_unlink(get_unique_path('tab')) as mutation_table_path:
        print_mutation_dict_to_file(td, mutation_table_path)
        plate_data_dict = get_cells_data_dict(srs, group_of_cell=group_of_cell)
        with relaxed_unlink(get_unique_path('tab')) as cell_data_path:
            write_cell_data_dict_to_file(plate_data_dict, cell_data_path)
            yield td, mutation_table_path, cell_data_path


@contextlib.contextmanager
def prep_plot_mutation_map_and_cell_data(srs, full_td, pl1, pl2, group_of_cell=lambda cell: cell.name):
    """
    The function
    :param srs: List of SampleReads
    :param full_td: Dict of cells with their loci length
    :param pl1:
    :param pl2:
    :param group_of_cell: list of the cell groups
    :return: (yield) td -
            mutation_table_path -
            cell_data_path -
            plot -
            output_plot -
    """
    with prep_mutation_map_and_cell_data(srs, full_td, pl1, pl2, group_of_cell=group_of_cell) \
            as (td, mutation_table_path, cell_data_path):
        with relaxed_unlink(get_unique_path('png')) as output_plot:
            plot = partial(
                basic_tree_enrichment_and_plotting,
                cell_data_path=cell_data_path,
                mutation_table_path=mutation_table_path,
                output_plot=output_plot,
                )
        yield td, mutation_table_path, cell_data_path, plot, output_plot
