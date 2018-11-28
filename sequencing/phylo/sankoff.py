from sequencing.calling.queries.mutation_maps import transpose_dict
import numpy as np
import pandas as pd

def get_probs(path='/home/dcsoft/s/ron/transition_table_exp9_ac_exvivo.csv'):
    names_lst = [i for i in range(1,31)]
    # prob_df = pd.read_csv('/home/dcsoft/s/ron/transition_table_nomodel_ac_exvivo.csv',names=names_lst )
    prob_df = pd.read_csv(path,names=names_lst )
    prob_df.index = np.arange(1, len(prob_df) + 1)
    return prob_df


def prep_tree(tree):
    for node in tree.leaf_node_iter():
        if node.taxon is not None:
            node.label = node.taxon.label
    try:
        tree.prune_taxa_with_labels(["root"])
        tree.taxon_namespace.remove_taxon_label('root')
    except KeyError:
        print("No taxon label root in tree")
    except LookupError:
        print("No taxon label root in tree")


def sankoff_node(node, nodes_probe_dict, prob_dict):
    if node.is_leaf():
        return
    child_nodes = list(node.child_nodes())
    for child_node in child_nodes:
        sankoff_node(child_node, nodes_probe_dict, prob_dict)
    node_probe_dict = sankoff_join(child_nodes, prob_dict, nodes_probe_dict)
    if node_probe_dict:
        nodes_probe_dict[node.label] = node_probe_dict


def sankoff_join(nodes, prob_dict, mutations_dict_of_single_loc):
    """Join two nodes together using the global cost matrix."""
    this = {}
    if any(node.label not in mutations_dict_of_single_loc for node in nodes):
        return this

    for i in prob_dict:  # possible genotypes
        #         min_left = np.inf
        #         min_right = np.inf
        min_nodes = dict()
        for node in nodes:
            min_node = np.inf
            for j in mutations_dict_of_single_loc[node.label]:
                this_cost = prob_dict[i][j] + mutations_dict_of_single_loc[node.label][j]
                min_node = min(min_node, this_cost)
            min_nodes[node] = min_node
        #         for j in mutations_dict_of_single_loc[node1.label]:
        #             this_cost = prob_dict[i][j] + mutations_dict_of_single_loc[node1.label][j]
        #             min_left = min(min_left, this_cost)

        #         for j in mutations_dict_of_single_loc[node2.label]:
        #             this_cost = prob_dict[i][j] + mutations_dict_of_single_loc[node2.label][j]
        #             min_right = min(min_right, this_cost)

        #         if np.isinf(min_left + min_right):
        #             continue
        if any(np.isinf(min_node) for min_node in min_nodes.values()):
            continue
        #         this[i] = min_left + min_right
        this[i] = sum(min_nodes.values())
    return this


def infer_loc_values_from_nodes_probe_dict(nodes_probe_dict):
    loc_values = dict()
    for n in nodes_probe_dict:
        if len(nodes_probe_dict[n]) == 1:
            loc_values[n] = list(nodes_probe_dict[n].keys())[0]
            continue
        if not nodes_probe_dict[n]:
            continue
        min_k, min_v = sorted(nodes_probe_dict[n].items(), key=lambda x: x[1])[0]
        loc_values[n] = min_k
    return loc_values


def calc_edge_len_noa_log_subtrees_no_divide(prob_dict, ud, vd,min_intersection_size=50):
    mutual_keys = ud.keys()&vd.keys()
    running_sum = 0
    if len(mutual_keys) < min_intersection_size:
        return None
    for loc in mutual_keys:
        running_sum += -np.log(prob_dict[ud[loc]][vd[loc]])
    return running_sum


def sankoffize_intermediate_nodes(
        tree, x_mutations_mono, transition_table='/home/dcsoft/s/ron/transition_table_exp9_ac_exvivo.csv'):
    prob_df = get_probs(path=transition_table)
    iprob_df = prob_df.applymap(lambda x: 1 - x)
    iprob_dict = iprob_df.to_dict()

    full_nodes_values_map = dict()

    for loc in x_mutations_mono:
        nodes_probe_dict = {k: {v: 0} for k, v in x_mutations_mono[loc].items()}
        sub_tree = tree.extract_tree_with_taxa_labels(nodes_probe_dict.keys())
        sankoff_node(sub_tree.seed_node, nodes_probe_dict, iprob_dict)  # updates nodes_probe_dict
        loc_values = infer_loc_values_from_nodes_probe_dict(nodes_probe_dict)
        #     infer_full_tree_from_subtree(tree, sub_tree, loc_values)  # updates loc_values
        full_nodes_values_map[loc] = loc_values

    # printing with different edge len
    trans_full_nodes_values_map = transpose_dict(full_nodes_values_map)
    return trans_full_nodes_values_map


def sankoffize_by_mono_loci(
        tree, x_mutations_mono, transition_table='/home/dcsoft/s/ron/transition_table_exp9_ac_exvivo.csv'):
    trans_full_nodes_values_map = sankoffize_intermediate_nodes(tree, x_mutations_mono)
    #     sankoffized_root = trans_full_nodes_values_map[None]
    func = calc_edge_len_noa_log_subtrees_no_divide
    # assign lengths
    prob_dict = get_probs(path=transition_table).to_dict()
    nodes_to_remove = []
    for e in tree.edges():
        if e.head_node is None or e.tail_node is None:
            continue
        if e.head_node == tree.seed_node or e.tail_node == tree.seed_node:
            continue
        if e.tail_node.label not in trans_full_nodes_values_map or e.head_node.label not in trans_full_nodes_values_map:
            dist = None
        else:
            dist = func(prob_dict, trans_full_nodes_values_map[e.tail_node.label],
                        trans_full_nodes_values_map[e.head_node.label], 15)
        if dist is None:
            nodes_to_remove.append(e.head_node)
            continue
        e.length = dist
    tree.ladderize(ascending=False)
    return tree


