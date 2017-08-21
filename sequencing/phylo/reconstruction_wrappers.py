import dendropy
from sequencing.phylo.plotting import prep_plot_mutation_map_and_cell_data
from misc.utils import relaxed_unlink, get_unique_path
from sequencing.phylo.utils import normalize_newick
from sequencing.phylo.triplets_wrapper import run_sagis_triplets
from sequencing.phylo.igor_wrapper import calculate_igors_tree
from sequencing.phylo.matlab_wrappers import distance_calculation, distance_based_reconstruction
from sequencing.phylo.tree_assessment.triplets_with_reference import triplets_scores_wrapper


def triplets_stats(srs, td, cells_group_map, pl1=50, pl2=50,
                   triplet_threshold=0.0, scoring_method='uri10', tripletsnumber=500000,
                   loci_filter='ncnr',
                   cells_to_be_used_as_root=tuple(['Ave']),
                   reference_tree_with_ids=None, tns_id_labels=None,
                   n=10000, min_d=3,
                  ):
    with prep_plot_mutation_map_and_cell_data(
        srs, td, pl1, pl2, group_of_cell=lambda cell: cells_group_map[cell.id]
    ) as (td, mutation_table_path, cell_data_path, plot):
        with relaxed_unlink(get_unique_path('newick')) as triplets_tree_path:
            run_sagis_triplets(
                textual_mutation_dict=td,
                cell_data_path=cell_data_path,
                cells_to_be_used_as_root=cells_to_be_used_as_root,
                newick_tree_path=triplets_tree_path,
                tripletsnumber=tripletsnumber,
                score_threshold=triplet_threshold,
                scoring_method=scoring_method,
                loci_filter=loci_filter)
            normalize_newick(triplets_tree_path, triplets_tree_path)
            triplets_tree = dendropy.Tree.get_from_path(
                triplets_tree_path,
                "newick",
                taxon_namespace=tns_id_labels)
            triplets_tree.encode_bipartitions()
            clustering_dict = plot(triplets_tree_path)
            if reference_tree_with_ids is not None:
                triplets_dict = triplets_scores_wrapper(
                    triplets_tree_path, reference_tree_with_ids, tns_id_labels, n=n, min_d=min_d)
                return clustering_dict, triplets_dict
            return clustering_dict


def igor_stats(srs, td, cells_group_map, pl1=50, pl2=50,
               reference_tree_with_ids=None, tns_id_labels=None,
               n=10000, min_d=3):
    with prep_plot_mutation_map_and_cell_data(
            srs, td, pl1, pl2,
            group_of_cell=lambda cell: cells_group_map[cell.id]) as (td, mutation_table_path, cell_data_path, plot):
        with relaxed_unlink(get_unique_path('newick')) as igor_tree_path:
            calculate_igors_tree(mutation_table_path, igor_tree_path, dist="ig1")
            normalize_newick(igor_tree_path, igor_tree_path)
            igor_tree = dendropy.Tree.get_from_path(
                igor_tree_path,
                "newick",
                taxon_namespace=tns_id_labels)
            res = igor_tree.encode_bipartitions()
            clustering_dict = plot(igor_tree_path)
            if reference_tree_with_ids is not None:
                triplets_dict = triplets_scores_wrapper(
                    igor_tree_path, reference_tree_with_ids, tns_id_labels, n=n, min_d=min_d)
                return clustering_dict, triplets_dict
            return clustering_dict


def nj_stats(srs, td, cells_group_map, pl1=50, pl2=50, distance_metric='ABS',
             reference_tree_with_ids=None, tns_id_labels=None, n=10000, min_d=3):
    with prep_plot_mutation_map_and_cell_data(
        srs, td, pl1, pl2, group_of_cell=lambda cell: cells_group_map[cell.id]) as (td, mutation_table_path, cell_data_path, plot):
        with relaxed_unlink(get_unique_path('tab')) as distance_matrix_path:
            distance_calculation(mutation_table_path, distance_matrix_path, distance_metric=distance_metric)
            with relaxed_unlink(get_unique_path('newick')) as nj_tree_path:
                with relaxed_unlink(get_unique_path('mat')) as tree_matlab_file:
                    distance_based_reconstruction(distance_matrix_path, nj_tree_path, tree_matlab_file)
                    normalize_newick(nj_tree_path, nj_tree_path)
                    nj_tree = dendropy.Tree.get_from_path(
                        nj_tree_path,
                        "newick",
                        taxon_namespace=tns_id_labels)
                    res = nj_tree.encode_bipartitions()
                    clustering_dict = plot(nj_tree_path)
                    if reference_tree_with_ids is not None:
                        triplets_dict = triplets_scores_wrapper(
                            nj_tree_path, reference_tree_with_ids, tns_id_labels, n=n, min_d=min_d)
                        return clustering_dict, triplets_dict
                    return clustering_dict
