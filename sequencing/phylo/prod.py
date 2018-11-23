import sys
sys.path.append("/home/dcsoft/clineage/")
import clineage.wsgi
from sequencing.calling.models import HighestPeaksProximityRatioFilteredBiSimCorSchemeModel, HighestPeaksMonoSimCorSchemeModel
from sequencing.runs.models import NGSRun, Demultiplexing
from sampling.models import Individual
from lib_prep.workflows.models import MagicalOM6BarcodedContent, MagicalPCR1BarcodedContent
import dendropy
from misc.utils import unlink, relaxed_unlink, get_unique_path
from sequencing.calling.queries.mutation_maps import get_root_genotypes
from sequencing.calling.queries.mutation_maps import transpose_dict
from lib_prep.workflows.models import MagicalOM6BarcodedContent
from sequencing.analysis.models import SampleReads, AdamHistogram, FullMSVHistogram
from sequencing.phylo.sankoff import prep_tree
import copy
from misc.utils import unlink, get_unique_path
from sequencing.phylo.triplets_wrapper import run_sagis_triplets_binary, convert_names_in_sagis_newick

from sequencing.phylo.utils import fix_directories, get_cells_group_map, textualize_d
import sys
sys.path.append('/home/dcsoft/clineage-simulation/')
from reconstruct import simplified_triplets_calculation



def get_bcs_and_histogram_type(ind, run_name, aa_runs=('nsr4', 'nsr5')):
    if run_name in aa_runs:
        bcs = MagicalPCR1BarcodedContent.objects.filter(content__cell__individual=ind)
        histogram_class = AdamHistogram
    else:
        bcs = MagicalOM6BarcodedContent.objects.filter(content__cell__individual=ind)
        histogram_class = FullMSVHistogram
    return bcs, histogram_class


aa_runs=frozenset({'nsr4', 'nsr5'})  # generally constant unless future aa runs are added

def get_multi_run_srs(ind, run_names):
    multi_run_srs = []
    runs_only = {run_name for run_name, demux_scheme_name in run_names}
    assert runs_only & aa_runs != {} or runs_only - aa_runs == runs_only # same individual can't have runs analyzed differently
    for run_name, demux_scheme_name in run_names:
        run = NGSRun.objects.get(name=run_name)
        demux = Demultiplexing.objects.get(ngs_run=run, demux_scheme__name=demux_scheme_name)
        bcs, histogram_class = get_bcs_and_histogram_type(ind, run_name, aa_runs=aa_runs)
        srs = SampleReads.objects.filter(demux=demux).filter(barcoded_content__in=bcs)
#         srs = [sr for sr in srs if sr.cell.composition == sc]
        skiped = [sr for sr in srs if not sr.cell.custom_group_labeling]
        if skiped:
            print('WARNING: removing {} cells with no group label'.format(len(skiped)))
        srs = [sr for sr in srs if sr.cell.custom_group_labeling]
        if demux_scheme_name == 'demux_scheme':
            srs = [sr for sr in srs if sr.barcoded_content.subclass.amplified_content.cell.custom_group_labeling in
                   ['Kera','Bulk']
                  ]
        multi_run_srs += srs
    return multi_run_srs, histogram_class


import numpy as np
def filter_mutation_map(mutations_dict, pl1=50, pl2=50):
    l1_counts = {kl1: len(mutations_dict[kl1]) for kl1 in mutations_dict}
    trans_mutations_dict = transpose_dict(mutations_dict)
    l2_counts = {kl2: len(trans_mutations_dict[kl2]) for kl2 in trans_mutations_dict}
    l1_lim = np.percentile(list(l1_counts.values()), pl1)
    l2_lim = np.percentile(list(l2_counts.values()), pl2)
    filtered_by_l1 = {kl1: d for kl1, d in mutations_dict.items() if l1_counts[kl1] >= l1_lim}
    filtered_by_l1_l2 = transpose_dict({kl2: d for kl2, d in transpose_dict(filtered_by_l1).items() if l2_counts[kl2] >= l2_lim})
    return filtered_by_l1_l2


def prep_mutation_table(multi_run_srs, d_by_ms, order=None, p_mono=0.7, filtering_percentiles=None):
    cells_group_map = get_cells_group_map(multi_run_srs)
    if filtering_percentiles is not None:
        d_by_ms = filter_mutation_map(d_by_ms, *filtering_percentiles)
    root = get_root_genotypes(d_by_ms, multi_run_srs, order=order, p_mono=p_mono)
    d_with_root = transpose_dict(copy.deepcopy(d_by_ms))
    d_with_root['root'] = root
    textual_d = textualize_d(d_with_root)
    return transpose_dict(textual_d), cells_group_map


def get_root_approximation_order(multi_run_srs, custom_group_labels):
    order = []
    for custom_group_labeling in custom_group_labels:
        srs = [sr for sr in multi_run_srs if sr.barcoded_content.subclass.amplified_content.cell.custom_group_labeling == custom_group_labeling]
        order.append(srs)
    return order


def sagi_wrapper(path_triplets_list_raw, triplets_tree_path, cell_id_map_for_sagi):
    # run sagis triplets binary
    # the output newick will have the numeric ids which correspond to the actual cell names
    with relaxed_unlink(get_unique_path('newick')) as index_labeled_output_newick:
        ret_code = run_sagis_triplets_binary(
            path_triplets_list_raw,
            index_labeled_output_newick,)

        # convert the numeric ids in newick to the actual cell names
        convert_names_in_sagis_newick(
            index_labeled_output_newick,
            triplets_tree_path,
            cell_id_map_for_sagi
        )


def end_to_end_triplets(
    textual_d_with_root_by_loc,
    output_prefix,
    choosing_method='mms',
    scoring_method='uri10',
    triplets_generator_name='splitable',
    # mono, splitable, full_bi, combined_likelihood, splitable_then_full, full_then_splittable
    loci_filter='ncnr',
    score_threshold=0,
):
    path_triplets_list_raw = output_prefix + '.triplets'
    cell_id_map_for_sagi = simplified_triplets_calculation(
        textual_d_with_root_by_loc,
        triplets_file=path_triplets_list_raw,
        triplets_generator_name=triplets_generator_name,
        score_threshold=score_threshold,
        choosing_method=choosing_method,
        scoring_method=scoring_method,
        printscores=True,
        loci_filter=loci_filter,
        tripletsnumber=50000000,
    )
    with unlink(get_unique_path('newick')) as triplets_tree_path:
        sagi_wrapper(path_triplets_list_raw, triplets_tree_path, cell_id_map_for_sagi)
        tree = dendropy.Tree.get_from_path(
            triplets_tree_path,
            "newick", )
    prep_tree(tree)

    normalized_newick_path = output_prefix + '.newick'
    tree.ladderize(ascending=False)
    tree.write_to_path(normalized_newick_path, 'newick', suppress_rooting=True)
    tree = dendropy.Tree.get_from_path(
        normalized_newick_path,
        "newick", )
    prep_tree(tree)