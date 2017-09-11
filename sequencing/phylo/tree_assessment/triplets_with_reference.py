import itertools
from sequencing.phylo.tree_assessment.utils import memory_expensive_random_choose
from dendropy import SeedNodeDeletionException
from frogress import bar
import dendropy


def challenge_triplets_generator(ref_tree, n=1000, min_d=3, exclude_direct_siblings=False):
    nodes = [n.taxon for n in ref_tree.leaf_nodes()]
    ndm = ref_tree.node_distance_matrix()
    pdm = ref_tree.phylogenetic_distance_matrix()
    for nodes_triplet in memory_expensive_random_choose(nodes, 3, n=n):
        pairs_in_triplet = list(itertools.combinations(nodes_triplet, 2))
        ref_close_pair = set(sorted(pairs_in_triplet, key=lambda x: pdm.path_edge_count(x[0],x[1]))[0])
        d = ndm.distance(ref_tree.mrca(taxa=nodes_triplet), ref_tree.mrca(taxa=ref_close_pair))
        if d < min_d:
            continue

        if exclude_direct_siblings:
            taxon_a, taxon_b = ref_close_pair
            distance_a_to_ab_mrca = ndm.distance(ref_tree.find_node_for_taxon(taxon_a), ref_tree.mrca(taxa=ref_close_pair), is_weighted_edge_distances=False)  # path_edge_count
            distance_b_to_ab_mrca = ndm.distance(ref_tree.find_node_for_taxon(taxon_b), ref_tree.mrca(taxa=ref_close_pair), is_weighted_edge_distances=False)
            if min(distance_a_to_ab_mrca, distance_b_to_ab_mrca) < 2:  # Exclude direct brothers
                continue
        yield d, nodes_triplet, ref_close_pair


def check_triplet(tree, nodes_triplet, ref_close_pair):
    try:
        rec_triplet = tree.extract_tree_with_taxa(nodes_triplet)
    except SeedNodeDeletionException:
        return None
    if len(rec_triplet.leaf_nodes()) < 3:
        return None
    rec_pdm = tree.phylogenetic_distance_matrix()
    pairs_in_triplet = list(itertools.combinations(nodes_triplet, 2))
    rec_close_pair = set(sorted(pairs_in_triplet, key=lambda x: rec_pdm.path_edge_count(x[0], x[1]))[0])  # you get the closest
    return ref_close_pair == rec_close_pair


def triplets_score(tree_in_question, reference_tree, n=1000, min_d=3):
    results_by_difficulty = dict()
    for d, nodes_triplet, ref_close_pair in bar(challenge_triplets_generator(reference_tree, n=n, min_d=min_d)):
        check_status = check_triplet(tree_in_question, nodes_triplet, ref_close_pair)
        if d not in results_by_difficulty:
            results_by_difficulty[d] = {True: 0, False: 0, None: 0}
        results_by_difficulty[d][check_status] += 1
    return results_by_difficulty


def triplets_scores_wrapper(rec_tree_path, reference_tree_with_ids, tns_id_labels, n=10000, min_d=3):
    rec_tree = dendropy.Tree.get_from_path(
        rec_tree_path,
        "newick",
        taxon_namespace=tns_id_labels)
    res = rec_tree.encode_bipartitions()
    return {
        k: v[True]/(v[True]+v[False]) if v[True]+v[False]>0 else None for k, v in triplets_score(
        rec_tree, reference_tree_with_ids, n=n, min_d=min_d).items()
    }
