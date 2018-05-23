from sequencing.phylo.matlab_wrappers import add_calculated_root_to_mutation_matrix
from sequencing.phylo.csv_writers import print_mutation_dict_to_file
from misc.utils import unlink, relaxed_unlink, get_unique_path
from sequencing.calling.queries.mutation_maps import transpose_dict
import os
import csv
import networkx as nx
import numpy as np
from frogress import bar
from plumbum import local
import dendropy
import sys
sys.path.append('/home/dcsoft/s/Ofir/triplets/triplets/')
from TMC_CLI import parse_mutations_table, paired_triplets_generator, format_triplet


def add_root_to_dict(
        textual_mutation_dict,
        cells_to_be_used_as_root,
):
    new_d = dict()
    if cells_to_be_used_as_root == ['Ave']:
        cells_to_be_used_as_root = set(c for loc in textual_mutation_dict for c in textual_mutation_dict[loc])
    else:
        cells_to_be_used_as_root = set(cells_to_be_used_as_root)
    root_collection = dict()
    for loc in textual_mutation_dict:
        for c in textual_mutation_dict[loc].keys() & cells_to_be_used_as_root:
            root_collection.setdefault(loc, []).append(textual_mutation_dict[loc][c])
    for loc in root_collection:
        val = int(np.median(root_collection[loc]))
        new_d.setdefault(loc, dict())['root'] = val
        new_d[loc].update(textual_mutation_dict[loc])
    return new_d


def map_cell_ids_for_sagi(rtd):
    rtd_for_sagi = dict()
    cell_id_map_for_sagi = dict()
    for i, cell_id in enumerate(rtd.keys()):
        cell_id_map_for_sagi[cell_id] = i
        if 'root' in cell_id:
            rtd_for_sagi[cell_id] = rtd[cell_id]
            continue
        rtd_for_sagi[i] = rtd[cell_id]
    return rtd_for_sagi, cell_id_map_for_sagi


def get_cells_and_root(mutation_table_path_for_triplets):
    calling = parse_mutations_table(mutation_table_path_for_triplets, inverse=True)
    # Verify the presence of a root cell in the input data.
    possible_roots = [cell for cell in calling if 'root' in cell]
    assert len(possible_roots) == 1
    root_label = possible_roots[0]
    root = (root_label, calling[root_label])
    cells = []
    for cell in calling:
        if cell == root_label:
            continue
        cells.append((cell, calling[cell]))
    assert len(cells) > 2
    return root, cells


def calculate_triplets_tree(
        textual_mutation_dict,
        triplets_file,
        cells_to_be_used_as_root=['Ave'],
        score_threshold=0,  # print scores
        choosing_method='mms',
        scoring_method='uri10',
        printscores=True,
        loci_filter='ncnr',
        sabc=0,
        tripletsnumber=5000,
):
    rtd = transpose_dict(textual_mutation_dict)
    if 'root' not in rtd:
        rtd = add_root_to_dict(
            textual_mutation_dict=textual_mutation_dict,
            cells_to_be_used_as_root=cells_to_be_used_as_root)
        rtd = transpose_dict(rtd)
    rtd_for_sagi, cell_id_map_for_sagi = map_cell_ids_for_sagi(rtd)
    with relaxed_unlink(get_unique_path("tab")) as mutation_table_path_for_triplets:
        print_mutation_dict_to_file(rtd_for_sagi, mutation_table_path_for_triplets)  # This file operation is redundent
        root, cells = get_cells_and_root(mutation_table_path_for_triplets)
    G = nx.Graph()
    G.add_nodes_from(cells)
    with open(triplets_file, 'w') as f:
        for triplet, pair, score in bar(paired_triplets_generator(root,
                                                              G,
                                                              loci_filter=loci_filter,
                                                              scoring_method=scoring_method,
                                                              choosing_method=choosing_method,
                                                              threshold=score_threshold,
                                                              triplets_num=tripletsnumber,
                                                              sabc=sabc)):
            f.write(format_triplet(triplet, pair, score, print_scores=printscores))
    return cell_id_map_for_sagi


def convert_tree_with_cell_id_map(taxon_name_space, cell_id_map_for_sagi):
    reverse_cell_id_map_for_sagi = {v: k for k, v in cell_id_map_for_sagi.items()}
    label_taxon_map = taxon_name_space.label_taxon_map()
    for l in label_taxon_map:
        original_label = reverse_cell_id_map_for_sagi[int(l)]
        taxa = label_taxon_map[l]
        taxa.label = original_label


def convert_names_in_sagis_newick(indexes_as_labels_tree, labeled_tree, cell_id_map_for_sagi):
    tns_sagi = dendropy.TaxonNamespace()
    indexes_tree = dendropy.Tree.get_from_path(
        indexes_as_labels_tree,
        "newick",
        taxon_namespace=tns_sagi)
    convert_tree_with_cell_id_map(tns_sagi, cell_id_map_for_sagi)
    indexes_tree.write_to_path(labeled_tree, 'newick', suppress_rooting=True)


tmc = local["/home/dcsoft/s/Noa/Tree_Analysis_2016/TMC/treeFromTriplets"]


class EmptyTripletsFile(Exception):
    pass


def run_sagis_triplets_binary(triplets_file, output_newick):
    # tmc(
    #     "-fid", triplets_file,
    #     "-frtN", output_newick)
    if os.stat(triplets_file).st_size == 0:
        raise EmptyTripletsFile('Empty file: '.format(triplets_file))
    sagi_cli = '/home/dcsoft/s/Noa/Tree_Analysis_2016/TMC/treeFromTriplets -fid {} -frtN {} -w 1 -index 2'.format(
        triplets_file, output_newick)
    return os.system(sagi_cli)


def run_sagis_triplets(
    textual_mutation_dict,
    cells_to_be_used_as_root,
    newick_tree_path,
    score_threshold=0,  # print scores
    choosing_method='mms',
    scoring_method='uri10',
    printscores=True,
    loci_filter='ncnr',
    sabc=0,
    tripletsnumber=5000,
):
    with relaxed_unlink(get_unique_path('triplets')) as triplets_list_path:
        cell_id_map_for_sagi = calculate_triplets_tree(
            textual_mutation_dict,
            triplets_list_path,
            cells_to_be_used_as_root,
            score_threshold=score_threshold,  # print scores
            choosing_method=choosing_method,
            scoring_method=scoring_method,
            printscores=printscores,
            loci_filter=loci_filter,
            sabc=sabc,
            tripletsnumber=tripletsnumber,
        )
        with relaxed_unlink(get_unique_path('newick')) as index_labeled_output_newick:
            ret_code = run_sagis_triplets_binary(triplets_list_path, index_labeled_output_newick)
            convert_names_in_sagis_newick(index_labeled_output_newick, newick_tree_path, cell_id_map_for_sagi)
