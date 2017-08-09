import csv
from misc.utils import unlink, get_unique_path
import matlab.engine
eng = matlab.engine.start_matlab()
eng.addpath(r'/home/dcsoft/s/Ofir/noa_matlab/Code/')
eng.addpath(r'/home/dcsoft/s/Ofir/noa_matlab/Code/General_Functions/')
eng.addpath(r'/home/dcsoft/s/Ofir/noa_matlab/Code/General_Functions/jsonlab-1.2/jsonlab/')


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
    eng.Distance_Calculation(
        'MStableInputFile', mutation_table_path,  # Input file (mutations table)
        'DistanceMatOutputFile', distance_matrix_output_path,  # Output file (Distance table)
        'DistMetric', distance_metric_map[distance_metric],  # METRIC_DICT = struct('ABS',0,'NORM_ABS',1,'EUC',2,'SQUARED',3,'MU_COUNT',4,'ML',5,'IGOR',6);
        'UseDataToReconstructTree', data_for_reconstruction,  # 'MS' / 'SNP' / 'MS_SNP'
        'CellsToBeAnalysed', filter_cells_by,  # Filter cells with data in 'MS' / 'SNP' / 'MS_SNP'  !!!Use Only MS until understanding this
        'MSweight', ms_snp_weight  # Relevant only when UseDataToReconstructTree='MS_SNP'
    )


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
    eng.calculate_root(
        cell_data_path,
        mutation_table_path,
        SNP_TAB_FILE,
        added_root_mutation_table_path,
        SNP_tab_output_file_with_Root,
        cells_to_be_used_as_root,
        JSON_duplicates_file_name
    )


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
    eng.Tree_Distance_Based_Reconstruction(
        'DistanceMatInputFile', distance_matrix_path,
        'Algorithm', algorithm,
        'TreeOutputFileNewick', tree_newick_output,
        'TreeOutputFileMatlab', tree_matlab_file,
    )


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
):
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
    )


def basic_tree_enrichment_and_plotting(
        tree_newick,
        cell_data_path,
        mutation_table_path,
):
    with unlink(get_unique_path('json')) as json1:
        with unlink(get_unique_path('json')) as json2:
            with unlink(get_unique_path('json')) as json3:
                with unlink(get_unique_path('json')) as json4:
                    with unlink(get_unique_path('json')) as json5:
                        with unlink(get_unique_path('tab')) as clustering_metrics_file:
                            tree_enrichment_and_plotting(
                                tree_newick,
                                cell_data_path,
                                mutation_table_path,
                                snp_table_path='',
                                Lineage_Output_File=clustering_metrics_file,
                                JSON_duplicates_file_name=json1,
                                JSON_cluster_width_file_name=json2,
                                JSON_cluster_color_file_name=json3,
                                JSON_leaf_color_file_name=json4,
                                JSON_legend_file_name=json5,
                            )
                            clustering_scores = read_statistics_file(clustering_metrics_file)
    return clustering_scores
