import matlab.engine
eng = matlab.engine.start_matlab()
eng.addpath(r'/home/dcsoft/s/Ofir/noa_matlab/Code/')
eng.addpath(r'/home/dcsoft/s/Ofir/noa_matlab/Code/General_Functions/')
eng.addpath(r'/home/dcsoft/s/Ofir/noa_matlab/Code/General_Functions/jsonlab-1.2/jsonlab/');


def distance_calculation(
        mutation_table_path,
        distance_matrix_output_path,
        distance_metric=1,
        data_for_reconstruction='MS',
        filter_cells_by='MS',
        ms_snp_weight=0,
):
    """
    MATLAB wrapper for MStableInputFile
    Args:
        mutation_table_path: Input file (mutations table)
        distance_matrix_output_path: Output file (Distance table)
        distance_metric: METRIC_DICT = struct('ABS',0,'NORM_ABS',1,'EUC',2,'SQUARED',3,'MU_COUNT',4,'ML',5,'IGOR',6)
        data_for_reconstruction: 'MS' / 'SNP' / 'MS_SNP'
        filter_cells_by: Filter cells with data in 'MS' / 'SNP' / 'MS_SNP'
        ms_snp_weight: Relevant only when UseDataToReconstructTree='MS_SNP'
    """
    assert data_for_reconstruction in ['MS', 'SNP', 'MS_SNP']
    assert filter_cells_by in ['MS', 'SNP', 'MS_SNP']
    eng.Distance_Calculation(
        'MStableInputFile', mutation_table_path,  # Input file (mutations table)
        'DistanceMatOutputFile', distance_matrix_output_path,  # Output file (Distance table)
        'DistMetric', distance_metric,  # METRIC_DICT = struct('ABS',0,'NORM_ABS',1,'EUC',2,'SQUARED',3,'MU_COUNT',4,'ML',5,'IGOR',6);
        'UseDataToReconstructTree', data_for_reconstruction,  # 'MS' / 'SNP' / 'MS_SNP'
        'CellsToBeAnalysed', filter_cells_by,  # Filter cells with data in 'MS' / 'SNP' / 'MS_SNP'  !!!Use Only MS until understanding this
        'MSweight', ms_snp_weight  # Relevant only when UseDataToReconstructTree='MS_SNP'
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
