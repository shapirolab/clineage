import os
import csv
from misc.utils import unlink, relaxed_unlink, get_unique_path
from clineage.settings import NOA_MATLAB
from .ete3_draw_tree import tree_draw
from .utils import plot_display_tree

import matlab.engine


def start_matlab_eng():
    eng = matlab.engine.start_matlab()
    matlab_general_functions = os.path.join(NOA_MATLAB, r'General_Functions/')
    matlab_json = os.path.join(matlab_general_functions, r'jsonlab-1.2/jsonlab/')
    eng.addpath(NOA_MATLAB)
    eng.addpath(matlab_general_functions)
    eng.addpath(matlab_json)
    return eng


def stop_matlab_eng(eng):
    eng.quit()


def sr_label_func(sr):
    return 'ID{}'.format(sr)


def ms_label_func(ms):
    return 'LOC_{}'.format(ms)


def distance_calculation(
        mutation_table_path,
        distance_matrix_output_path,
        distance_metric='MU_COUNT',
        data_for_reconstruction='MS',
        filter_cells_by='MS',
        ms_snp_weight=0,
):
    """
    MATLAB wrapper for MStableInputFile
    Args:
        mutation_table_path: Input file (mutations table)
        distance_matrix_output_path: Output file (Distance table)
        distance_metric: METRIC_DICT = ['ABS','NORM_ABS','EUC','SQUARED','MU_COUNT','ML','IGOR']
        data_for_reconstruction: 'MS' / 'SNP' / 'MS_SNP'
        filter_cells_by: Filter cells with data in 'MS' / 'SNP' / 'MS_SNP'
        ms_snp_weight: Relevant only when UseDataToReconstructTree='MS_SNP'
    """
    assert data_for_reconstruction in ['MS', 'SNP', 'MS_SNP']
    assert filter_cells_by in ['MS', 'SNP', 'MS_SNP']
    distance_metric_map = {
        'ABS': 0,
        'NORM_ABS': 1,
        'EUC': 2,
        'SQUARED': 3,
        'MU_COUNT': 4,
        'ML': 5,
        'IGOR': 6}
    assert distance_metric in distance_metric_map.keys()
    eng = start_matlab_eng()
    eng.Distance_Calculation(
        'MStableInputFile', mutation_table_path,  # Input file (mutations table)
        'DistanceMatOutputFile', distance_matrix_output_path,  # Output file (Distance table)
        'DistMetric', distance_metric_map[distance_metric],  # METRIC_DICT = struct('ABS',0,'NORM_ABS',1,'EUC',2,'SQUARED',3,'MU_COUNT',4,'ML',5,'IGOR',6);
        'UseDataToReconstructTree', data_for_reconstruction,  # 'MS' / 'SNP' / 'MS_SNP'
        'CellsToBeAnalysed', filter_cells_by,  # Filter cells with data in 'MS' / 'SNP' / 'MS_SNP'  !!!Use Only MS until understanding this
        'MSweight', ms_snp_weight  # Relevant only when UseDataToReconstructTree='MS_SNP'
    )
    stop_matlab_eng(eng)


def add_calculated_root_to_mutation_matrix(
    mutation_table_path,
    cell_data_path,
    added_root_mutation_table_path,
    cells_to_be_used_as_root=['Ave'],
    SNP_TAB_FILE='',
    SNP_tab_output_file_with_Root='',
    JSON_duplicates_file_name='/dev/null',
):
    """
    MATLAB wrapper for calculate_root
    Args:
        mutation_table_path: 
        cell_data_path: 
        added_root_mutation_table_path: 
        cells_to_be_used_as_root: 
        SNP_TAB_FILE: 
        SNP_tab_output_file_with_Root: 
        JSON_duplicates_file_name: 

    Returns:

    """
    eng = start_matlab_eng()
    eng.calculate_root(
        cell_data_path,
        mutation_table_path,
        SNP_TAB_FILE,
        added_root_mutation_table_path,
        SNP_tab_output_file_with_Root,
        cells_to_be_used_as_root,
        JSON_duplicates_file_name
    )
    stop_matlab_eng(eng)


def distance_based_reconstruction(
        distance_matrix_path,
        tree_newick_output,
        tree_matlab_file='/dev/null',
        algorithm='NJ',
):
    """
    MATLAB wrapper for Tree_Distance_Based_Reconstruction
    Args:
        distance_matrix_path: 
        tree_newick_output: 
        tree_matlab_file: 
        algorithm: 

    Returns:

    """
    eng = start_matlab_eng()
    eng.Tree_Distance_Based_Reconstruction(
        'DistanceMatInputFile', distance_matrix_path,
        'Algorithm', algorithm,
        'TreeOutputFileNewick', tree_newick_output,
        'TreeOutputFileMatlab', tree_matlab_file,
    )
    stop_matlab_eng(eng)


def read_statistics_file(path):
    with open(path) as f:
        rdr = csv.DictReader(f, dialect='excel-tab')
        rows = list(rdr)
        assert len(rows) == 1
        row = rows[0]
        return row


def tree_enrichment_and_plotting(
        tree_newick,
        cell_data_path,
        mutation_table_path,
        snp_table_path='',
        Lineage_Output_File='/dev/null',
        JSON_duplicates_file_name='/dev/null',
        JSON_cluster_width_file_name='/dev/null',
        JSON_cluster_color_file_name='/dev/null',
        JSON_leaf_color_file_name='/dev/null',
        JSON_legend_file_name='/dev/null',
        JSON_bootstrap_file_name='/dev/null',
        JSON_leaf_order_file_name='/dev/null',
        JSON_leaf_label_file_name='/dev/null',
):
    eng = start_matlab_eng()
    eng.Draw_tree_with_enrichment(
        cell_data_path,  # CELL_DATA_FILE,
        mutation_table_path,  # MS_TAB_FILE,  # this is actually a list of the cells to filter by
        snp_table_path,  # SNP_TAB_FILE,  # color by snp values in table
        tree_newick,  # Tree_input_file_newick
        Lineage_Output_File,  # general statistics file (doesn't really work at all)
        JSON_duplicates_file_name,
        JSON_cluster_width_file_name,
        JSON_cluster_color_file_name,
        JSON_leaf_color_file_name,
        JSON_legend_file_name,
        JSON_bootstrap_file_name,
        JSON_leaf_order_file_name,
        JSON_leaf_label_file_name,
    )
    stop_matlab_eng(eng)


def basic_tree_enrichment_and_plotting(
        tree_newick,
        cell_data_path,
        mutation_table_path,
        output_plot,
        tree_name='Tree',
        tree_scale='linear',
        tree_rotation=True,
        font_size=7,
        font_legend=7,
        node_size=3,
        fig_width=500,
        fig_height=300,
        fig_dpi=120,
        scale_rate=None,
        distance_factor=1,
        y_scale=True,


):
    with relaxed_unlink(get_unique_path('json')) as duplicates_file:
        with relaxed_unlink(get_unique_path('json')) as cluster_width_file:
            with relaxed_unlink(get_unique_path('json')) as cluster_color_file:
                with relaxed_unlink(get_unique_path('json')) as leaf_color_file:
                    with relaxed_unlink(get_unique_path('json')) as legend_file:
                        with relaxed_unlink(get_unique_path('tab')) as clustering_metrics_file:
                            with relaxed_unlink(get_unique_path('json')) as bootstrap_file:
                                with relaxed_unlink(get_unique_path('json')) as leaf_order_file:
                                    with relaxed_unlink(get_unique_path('json')) as leaf_label_file:

                                        tree_enrichment_and_plotting(
                                            tree_newick,
                                            cell_data_path,
                                            mutation_table_path,
                                            snp_table_path='',
                                            Lineage_Output_File=clustering_metrics_file,
                                            JSON_duplicates_file_name=duplicates_file,
                                            JSON_cluster_width_file_name=cluster_width_file,
                                            JSON_cluster_color_file_name=cluster_color_file,
                                            JSON_leaf_color_file_name=leaf_color_file,
                                            JSON_legend_file_name=legend_file,
                                            JSON_bootstrap_file_name=bootstrap_file,
                                            JSON_leaf_order_file_name=leaf_order_file,
                                            JSON_leaf_label_file_name=leaf_label_file,
                                        )

                                        # tree, tstyle = tree_draw(
                                        #     tree_newick,
                                        #     tree_name=tree_name,
                                        #     duplicate_file=duplicates_file,
                                        #     clustering_sizes_file=cluster_width_file,
                                        #     clustering_colors_file=cluster_color_file,
                                        #     cell_colors_file=leaf_color_file,
                                        #     legend_file=legend_file,
                                        #     intermediate_node_sizes_file=bootstrap_file,
                                        #     intermediate_node_labels_file=None,
                                        #     order_vector_file=leaf_order_file,
                                        #     leaf_labels_file=leaf_label_file,
                                        #     tree_scale=tree_scale,
                                        #     tree_rotation=tree_rotation,
                                        #     font_size=font_size,
                                        #     font_legend=font_legend,
                                        #     node_size=node_size,
                                        #     scale_rate=scale_rate,
                                        #     distance_factor=distance_factor,
                                        #     y_scale=y_scale,
                                        # )

                                        # plot_display_tree(
                                        #     tree,
                                        #     tstyle,
                                        #     output_plot,
                                        #     fig_width=fig_width,
                                        #     fig_height=fig_height,
                                        #     fig_dpi=fig_dpi,
                                        # )

                                        clustering_scores = read_statistics_file(clustering_metrics_file)
    return clustering_scores

