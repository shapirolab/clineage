__author__ = 'veronika'

import primer3
import re
import argparse
import sys
sys.path.append('/home/barakor/ViennaRNA/lib/python3.5/site-packages/')
import RNA

from genomes.models import DNASlice
from targeted_enrichment.planning.models import Target

GLOBAL_ARGS = {
    'PRIMER_TASK': 'pick_detection_primers',
    'PRIMER_OPT_SIZE': 23,
    'PRIMER_PICK_INTERNAL_OLIGO': 1,
    'PRIMER_INTERNAL_MAX_SELF_END': 8,
    'PRIMER_MIN_SIZE': 20,
    'PRIMER_MAX_SIZE': 27,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 57.0,
    'PRIMER_MAX_TM': 63.0,
    'PRIMER_TM_FORMULA': 1,
    'PRIMER_MIN_GC': 20,
    'PRIMER_MAX_GC': 80,
    'PRIMER_MAX_POLY_X': 100,
    'PRIMER_INTERNAL_MAX_POLY_X': 100,
    'PRIMER_DNA_CONC': 50.0,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_SELF_ANY': 12,
    'PRIMER_MAX_SELF_END': 8,
    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    'PRIMER_PAIR_MAX_COMPL_END': 8,
    'PRIMER_PAIR_MAX_DIFF_TM': 3.0,
    'PRIMER_NUM_RETURN': 10000,
    'PRIMER_FILE_FLAG': 0,
    'P3_FILE_FLAG': 0,
    'PRIMER_PRODUCT_SIZE_RANGE': [[130, 300]],
}


def create_amplicon_for_primer3(slice, slice_margin):
    template_slice = DNASlice(chromosome=slice.chromosome, start_pos=slice.start_pos-slice_margin, end_pos=slice.end_pos+slice_margin)
    return template_slice


def primer3_design(target, **kwargs):

    margins = kwargs.get('margins', 0)
    slice_margin = kwargs.get('slice_margin', 1000)
    primer3_template_margins = create_amplicon_for_primer3(target.slice, slice_margin)

    GLOBAL_ARGS.update(kwargs)
    GLOBAL_ARGS['PRIMER_PRODUCT_SIZE_RANGE'] = [[kwargs.get('min_size', 130), kwargs.get('max_size', 300)]]

    primer_details = dict()
    primer_details['SEQUENCE_ID'] = target.id
    primer_details['SEQUENCE_TEMPLATE'] = primer3_template_margins.sequence.seq.decode('utf-8')
    slice_len = len(primer_details['SEQUENCE_TEMPLATE'])-2*slice_margin
    primer_details['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [slice_margin -
                                                             (GLOBAL_ARGS['PRIMER_PRODUCT_SIZE_RANGE'][0][1] -
                                                              (slice_len//2+margins)),
                                                             GLOBAL_ARGS['PRIMER_PRODUCT_SIZE_RANGE'][0][1] -
                                                             (slice_len//2+margins),
                                                             slice_margin+slice_len,
                                                             GLOBAL_ARGS['PRIMER_PRODUCT_SIZE_RANGE'][0][1] -
                                                             (slice_len//2+margins)]

    primer3_output = primer3.bindings.designPrimers(primer_details, GLOBAL_ARGS)
    primer3_output['SEQUENCE_ID'] = target.id
    primer3_output['TARGET_START'] = (slice_margin, target.slice.end_pos - target.slice.start_pos)

    return primer3_output


def parse_primers(primer3_output):
    target_primers = dict()

    number_of_primers = len(primer3_output) // 30
    for i in range(1, number_of_primers):
        target_primers[i] = {}
        target_primers[i]['LEFT'] = primer3_output['PRIMER_LEFT_{}_SEQUENCE'.format(i)]
        target_primers[i]['RIGHT'] = primer3_output['PRIMER_RIGHT_{}_SEQUENCE'.format(i)]
        target_primers[i]['LEFT_TM'] = primer3_output['PRIMER_LEFT_{}_TM'.format(i)]
        target_primers[i]['RIGHT_TM'] = primer3_output['PRIMER_RIGHT_{}_TM'.format(i)]
        target_primers[i]['LEFT_START'] = primer3_output['PRIMER_LEFT_{}'.format(i)]
        target_primers[i]['RIGHT_START'] = primer3_output['PRIMER_RIGHT_{}'.format(i)]
    target_primers['TARGET_START'] = primer3_output['TARGET_START']
    target_primers['id'] = primer3_output['SEQUENCE_ID']
    return target_primers


def filter_by_size(filter_list, target_primers, best_size=160, margin_size=250):

    for tar in filter_list:
        tar_size = target_primers[tar]['RIGHT_START'][0] - target_primers[tar]['LEFT_START'][0]
        if not best_size + margin_size >= tar_size >= best_size - margin_size:
            filter_list.remove(tar)

    return filter_list


def filter_by_delta_g(filter_list, target_primers, target, relative_target_pos, chromosome, delta_min=-70, delta_max=0):
    for tar in filter_list:
        template_slice = DNASlice(chromosome=chromosome,
                                  start_pos=target.slice.start_pos - (relative_target_pos -
                                                                      target_primers[tar]['LEFT_START'][0]),
                                  end_pos=target.slice.start_pos + (target_primers[tar]['RIGHT_START'][0] -
                                                                    relative_target_pos)
                                  )
        ter = template_slice.sequence.seq.decode('utf-8')
        dctpd = RNA.fold(ter)[1]

        if not delta_min < dctpd <= delta_max:
            filter_list.remove(tar)

    return filter_list


def sort_best_primers(primer3_output, **kwargs):
    """
    filtering primers based on different parameters that we choose
    :param primer3_output: primers dict for targets
    :param size_filter: primers dict for targets
    :param deltag_filter: primers dict for targets
    :param best_size: primers dict for targets
    :param margin_size: primers dict for targets
    :param delta_min: primers dict for targets
    :param delta_max: primers dict for targets
    :return: chosen_target_primers, discarded_targets
    """
    size_filter = kwargs.get('size_filter', True)
    deltag_filter = kwargs.get('deltag_filter', False)

    target_primers = parse_primers(primer3_output)
    discarded_targets = []
    chosen_target_primers = []

    # looking for the best length amplicon for PCR
    target = Target.objects.get(id=target_primers['id'])
    relative_target_pos = target_primers['TARGET_START'][0]
    chromosome = target.slice.chromosome

    filter_list = list(set(target_primers.keys())-{'id', 'TARGET_START'})
    filter_list.sort()

    if size_filter:
        filter_list = filter_by_size(filter_list, target_primers,
                                     best_size=kwargs.get('best_size', 150), margin_size=kwargs.get('margin_size', 250))

    # checking deltaG for amplicons
    if deltag_filter:
        filter_list = filter_by_delta_g(filter_list, target_primers, target, relative_target_pos, chromosome,
                                        delta_min=kwargs.get('delta_min', -70), delta_max=kwargs.get('delta_max', 0))

    if filter_list:
        chosen_target_primers = target_primers[filter_list[0]]
        chosen_target_primers['id'] = target_primers['id']
        chosen_target_primers['TARGET_START'] = target_primers['TARGET_START']

    if not chosen_target_primers:
        discarded_targets = target_primers['id']

    return chosen_target_primers, discarded_targets
