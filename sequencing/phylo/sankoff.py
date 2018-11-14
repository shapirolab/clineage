from frogress import bar
from sequencing.calling.queries.mutation_maps import transpose_dict
from sequencing.analysis.snps.models import SNPReads
from sequencing.analysis.snps.parse_snps import snp_calling

from sequencing.calling.queries.mutation_maps import filter_mutation_map
from lib_prep.workflows.models import MagicalOM6BarcodedContent
from targeted_enrichment.amplicons.models import Amplicon
import functools
from collections import Counter
from frogress import bar
from sequencing.analysis.models import SampleReads, AdamMergedReads, AdamReadsIndex, AdamMarginAssignment, \
    AdamAmpliconReads, AdamHistogram, HistogramEntryReads, AdamMSVariations, \
    MicrosatelliteHistogramGenotype, HistogramEntryReads, Histogram
from targeted_enrichment.planning.models import Microsatellite
from sequencing.calling.models import CalledAlleles
from django.contrib.auth.models import User
import pandas as pd
from sequencing.analysis.full_msv.models import FullMSVHistogram
from sequencing.calling.hist import Histogram as dHistogram
from sequencing.calling.simcor.hist_analysis import better_get_far_apart_highest_peaks
from sequencing.calling.queries.mutation_maps import filter_mutation_map

import numpy as np


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