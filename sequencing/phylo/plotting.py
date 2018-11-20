import copy
import contextlib
from sequencing.calling.queries.mutation_maps import filter_bipartition_loci, filter_mutation_map
from sequencing.calling.queries.formatters import get_cells_data_dict
from sequencing.phylo.csv_writers import write_cell_data_dict_to_file, print_mutation_dict_to_file
from functools import partial
from sequencing.phylo.matlab_wrappers import basic_tree_enrichment_and_plotting
from misc.utils import unlink, relaxed_unlink, get_unique_path


@contextlib.contextmanager
def prep_mutation_map_and_cell_data(srs, full_td, pl1, pl2, group_of_cell=lambda cell: cell.name, complet_color_map_override=None):
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
        plate_data_dict = get_cells_data_dict(srs, group_of_cell=group_of_cell, complet_color_map_override=complet_color_map_override)
        with relaxed_unlink(get_unique_path('tab')) as cell_data_path:
            write_cell_data_dict_to_file(plate_data_dict, cell_data_path)
            yield td, mutation_table_path, cell_data_path


@contextlib.contextmanager
def prep_plot_mutation_map_and_cell_data(srs, full_td, pl1, pl2, group_of_cell=lambda cell: cell.name, complet_color_map_override=None):
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
    with prep_mutation_map_and_cell_data(srs, full_td, pl1, pl2, group_of_cell=group_of_cell, complet_color_map_override=complet_color_map_override) \
            as (td, mutation_table_path, cell_data_path):
        with relaxed_unlink(get_unique_path('png')) as output_plot:
            plot = partial(
                basic_tree_enrichment_and_plotting,
                cell_data_path=cell_data_path,
                mutation_table_path=mutation_table_path,
                output_plot=output_plot,
                )
        yield td, mutation_table_path, cell_data_path, plot, output_plot


import csv

@contextlib.contextmanager
def minimal_mutations_table_context_wrapper(tree):
    with relaxed_unlink(get_unique_path('tab')) as tmp_file_fake_mutations_table:
        print_mutation_dict_to_file(
            {node.taxon.label: {'FAKE_LOC': 12} for node in tree.leaf_nodes()},
            tmp_file_fake_mutations_table)
        yield tmp_file_fake_mutations_table


@contextlib.contextmanager
def minimal_cell_data_context_wrapper(tree, groups_map):
    """
    cell data file is indexed by integers which correspond to "ID{integer}" leaf labels in the newick
    Args:
        tree: dendropy tree object, "ID{integer}" leaf labels
        groups_map: dictionary mapping the ids (without "ID" prefix) to textual group names

    Returns:

    """
    with relaxed_unlink(get_unique_path('tab')) as tmp_file_fake_cell_data_table:
        with open(tmp_file_fake_cell_data_table, 'w') as f:
            fieldnames = ['Cell_Group', 'Cell_Group_short', 'Cell_ID', 'Sample_Reads_ID', ]
            writer = csv.DictWriter(f, fieldnames=fieldnames, dialect='excel-tab')
            writer.writeheader()
            for node in tree.leaf_nodes():
                writer.writerow(
                    {
                        'Cell_Group': groups_map[int(node.taxon.label[2:])],
                        'Cell_Group_short': '',
                        'Cell_ID': int(node.taxon.label[2:]),  # must be a number
                        'Sample_Reads_ID': int(node.taxon.label[2:])  # must be a number
                    })
        yield tmp_file_fake_cell_data_table


@contextlib.contextmanager
def minimal_noa_context_wrapper(tree, groups_map):
    with minimal_mutations_table_context_wrapper(tree) as tmp_file_fake_mutations_table:
        with minimal_cell_data_context_wrapper(tree, groups_map) as tmp_file_fake_cell_data_table:
            yield tmp_file_fake_mutations_table, tmp_file_fake_cell_data_table


def minimal_noa_plot(
        tree,
        groups_map,
        output_prefix,
        title='',
        node_size=5,
        font_legend=10,
        fig_width=800,
        fig_height=800,
        scale_rate=50,
        tree_rotation=True):
    """
    taxon labels should be unique, used for both "Cell_ID" and "Sample_Reads_ID" columns
    grouping map is indexed by taxon labels
    column "Cell_Group_short" is required by matlab but not actually used
    output_prefix - full path + filename prefix for output newick, png, pdf
    kwargs node_size, font_legend, fig_width, fig_height,  scale_rate, tree_rotation are only relevant
    for the python plotter that is currently disabled (hardcoded)
    """
    newick_path = output_prefix + '.newick'
    output_plot = output_prefix + '.png'
    tree.write_to_path(newick_path, 'newick', suppress_rooting=True)
    with minimal_noa_context_wrapper(tree, groups_map) as (tmp_mt, tmp_cdt):
        clustering_dict = basic_tree_enrichment_and_plotting(
            newick_path,
            cell_data_path=tmp_cdt,
            mutation_table_path=tmp_mt,
            output_plot=output_plot,
            tree_name=title,
            node_size=node_size,
            font_legend=font_legend,
            fig_width=fig_width,
            fig_height=fig_height,
            scale_rate=scale_rate,
            tree_rotation=tree_rotation)
        return clustering_dict
