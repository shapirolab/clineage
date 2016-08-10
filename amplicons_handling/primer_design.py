__author__ = 'veronika'

import primer3
import re
import argparse
import pysam
from collections import defaultdict
from collections import Counter
from frogress import bar
from plumbum import local
# import sys
# sys.path.append('/home/barakor/ViennaRNA/lib/python3.5/site-packages/')
# import RNA

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
    'PRIMER_MIN_GC': 20.0,
    'PRIMER_MAX_GC': 80.0,
    'PRIMER_MAX_POLY_X': 100,
    'PRIMER_INTERNAL_MAX_POLY_X': 100,
    'PRIMER_DNA_CONC': 50.0,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_SELF_ANY': 12,
    'PRIMER_MAX_SELF_END': 8,
    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    'PRIMER_PAIR_MAX_COMPL_END': 8,
    'PRIMER_PAIR_MAX_DIFF_TM': 3,
    'PRIMER_NUM_RETURN': 1000,
    'PRIMER_FILE_FLAG': 0,
    'P3_FILE_FLAG': 0,
    'PRIMER_PRODUCT_SIZE_RANGE': [[130, 400]]
}


def create_amplicon_for_primer3(slice, slice_margin):
    template_slice = DNASlice(chromosome=slice.chromosome, start_pos=slice.start_pos-slice_margin, end_pos=slice.end_pos+slice_margin)
    return template_slice


def primer3_design(target, min_size=130, max_size=400, primer_num_rerun=10000, opt_tm=60.0, max_tm=63.0,
                   min_tm=57.0, delta_tm=3.0, margins=0, slice_margin=1000):
    primer3_template_margins = create_amplicon_for_primer3(target.slice, slice_margin)

    GLOBAL_ARGS['PRIMER_NUM_RETURN'] = primer_num_rerun
    GLOBAL_ARGS['PRIMER_OPT_TM'] = opt_tm
    GLOBAL_ARGS['PRIMER_MIN_TM'] = min_tm
    GLOBAL_ARGS['PRIMER_MAX_TM'] = max_tm
    GLOBAL_ARGS['PRIMER_PAIR_MAX_DIFF_TM'] = delta_tm
    GLOBAL_ARGS['PRIMER_PRODUCT_SIZE_RANGE'] = [[min_size, max_size]]

    primer_details = dict()
    primer_details['SEQUENCE_ID'] = target.name
    primer_details['SEQUENCE_TEMPLATE'] = primer3_template_margins.sequence.seq.decode('utf-8')
    slice_len = len(primer_details['SEQUENCE_TEMPLATE'])-2*slice_margin
    primer_details['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [slice_margin-(max_size-(slice_len//2+margins)),
                                                             max_size-(slice_len//2+margins),
                                                             slice_margin+slice_len,
                                                             max_size-(slice_len//2+margins)]

    primer3_output = primer3.bindings.designPrimers(primer_details, GLOBAL_ARGS)
    primer3_output['SEQUENCE_ID'] = target.id
    primer3_output['TARGET_START'] = (slice_margin, target.slice.end_pos - target.slice.start_pos)

    return primer3_output


def parse_primer3_output(output_name):
    fasta_output_file = open(output_name, 'rb')
    primer3_output = fasta_output_file.read()
    fasta_string = primer3_output.split('\n=\n')
    pid = re.compile('SEQUENCE_ID=(\w+)')
    pl = re.compile('PRIMER_LEFT_\d+_SEQUENCE=([actgACTG]+)')
    pr = re.compile('PRIMER_RIGHT_\d+_SEQUENCE=([actgACTG]+)')
    target_primers = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
    for seq in fasta_string[:-1]:
        seq_id = pid.findall(seq)[0]
        left_primers = pl.findall(seq)
        right_primers = pr.findall(seq)
        assert len(left_primers) == len(right_primers)
        for i, primer_left in enumerate(left_primers):
            primer_right = right_primers[i]
            target_primers[seq_id][i]['LEFT'] = primer_left
            target_primers[seq_id][i]['RIGHT'] = primer_right
    return target_primers


def bowtie2_design(input_fasta_file, output_file, bowtie2_index, output_name):
    target_primers = parse_primer3_output(output_name)
    out_string = ''
    for seq_id in list(target_primers.keys()):
        for primer_number in list(target_primers[seq_id]['LEFT'].keys()):
            out_string += '>PRIMER_LEFT_{}_{}\n{}\n'.format(primer_number, seq_id, target_primers[seq_id][primer_number]['LEFT'])
            out_string += '>PRIMER_RIGHT_{}_{}\n{}\n'.format(primer_number, seq_id, target_primers[seq_id][primer_number]['RIGHT'])
    primer_data_check = '{}.fa'.format(str(input_fasta_file))
    # print('writing primers fasta file {}'.format(primer_data_check))
    with open(primer_data_check, 'w+') as primers_output:
        primers_output.write(out_string)
    sam_file = '{}.sam'.format(str(output_file))
    print('running bowtie sam file output:{}'.format(sam_file))

    BOWTIE_2 = local[settings.PRIMER3_PATH]
    BOWTIE_2('-k 2 {} -f -U {} -S {}'.format(bowtie2_index, primer_data_check, sam_file))

    # s = '{} -k 2 {} -f -U {} -S {}'.format(settings.BOWTIE2_PATH, bowtie2_index, primer_data_check, sam_file)
    # os.system(s)
    return sam_file, primer_data_check, target_primers


def primer_count_from_sam_file(sam_file):
    samfile = pysam.Samfile(sam_file, "r")
    primer_keys = []
    primers_names = defaultdict(lambda: defaultdict(int))
    for sam_primer in samfile:
        primer_keys.append(sam_primer.qname)
        primers_names[sam_primer.qname]['Start_pos'] = sam_primer.reference_start
        primers_names[sam_primer.qname]['End_pos'] = sam_primer.reference_end-1
    name_count = Counter(primer_keys)
    return name_count, primers_names


def sort_best_primers(sam_file, target_primers, margins=300):
    name_count, primers_names = primer_count_from_sam_file(sam_file)
    chosen_target_primers = defaultdict(lambda: defaultdict(str))
    discarded_targets = []
    for target_id in bar(list(target_primers.keys())):
        for primer_number in sorted(target_primers[target_id].keys()):
            bowtie2_key_left = 'PRIMER_LEFT_{}_{}'.format(primer_number['LEFT'], target_id)
            bowtie2_key_right = 'PRIMER_RIGHT_{}_{}'.format(primer_number['RIGHT'], target_id)
            right_primer = target_primers[target_id][primer_number]['RIGHT']
            left_primer = target_primers[target_id][primer_number]['LEFT']
            left_pos = primers_names[bowtie2_key_left]['Start_pos']
            right_pos = primers_names[bowtie2_key_right]['End_pos']
            if right_pos - left_pos <= margins:
                break

        if left_primer and right_primer:
            target = Target.objects.get(pk=target_id)
            # try:
            #     primer_fw, primer_rv = check_primers(target,
            #                                          left_primer,
            #                                          right_primer,
            #                                          target_enrichment_type=TargetEnrichmentType.objects.get(
            #                                              name='PCR_with_tails'),
            #                                          margins=margins)
            #     #test_new?
            #     chosen_target_primers[target_id]['LEFT'] = left_primer
            #     chosen_target_primers[target_id]['RIGHT'] = right_primer
            # except AmpliconCollisionError:
            #     discarded_targets.append(target_id)
            # except PrimerLocationError:
            #     discarded_targets.append(target_id)
            chosen_target_primers[target_id]['LEFT'] = left_primer
            chosen_target_primers[target_id]['RIGHT'] = right_primer
        else:
            discarded_targets.append(target_id)

    return chosen_target_primers, discarded_targets


def sort_unique_primers(sam_file, target_primers, margins=300):
    name_count, primers_names = primer_count_from_sam_file(sam_file)
    chosen_target_primers = defaultdict(lambda: defaultdict(str))
    discarded_targets = []
    for target_id in bar(list(target_primers.keys())):
        left_primer = None
        for primer_number_fr in sorted(target_primers[target_id]['LEFT'].keys()):
            bowtie2_key_left = 'PRIMER_LEFT_{}_{}'.format(primer_number_fr, target_id)
            if name_count[bowtie2_key_left] == 1:
                left_primer = target_primers[target_id][primer_number_fr]['LEFT']
                right_primer = None
                for primer_number_rv in sorted(target_primers[target_id]['RIGHT'].keys()):
                    bowtie2_key_right = 'PRIMER_RIGHT_{}_{}'.format(primer_number_rv, target_id)
                    if name_count[bowtie2_key_right] == 1:
                        right_primer = target_primers[target_id][primer_number_rv]['RIGHT']
                        left_pos = primers_names[bowtie2_key_left]['Start_pos']
                        right_pos = primers_names[bowtie2_key_right]['End_pos']
                        if right_pos - left_pos <= margins:
                            break
                        else:
                            right_primer = None

        if left_primer and right_primer:
            target = Target.objects.get(pk=target_id)
            try:
                primer_fw, primer_rv = check_primers(target,
                                 left_primer,
                                 right_primer,
                                 target_enrichment_type=TargetEnrichmentType.objects.get(name='PCR_with_tails'),
                                 margins=margins)
                chosen_target_primers[target_id]['LEFT'] = left_primer
                chosen_target_primers[target_id]['RIGHT'] = right_primer
            except AmpliconCollisionError:
                discarded_targets.append(target_id)
            except PrimerLocationError:
                discarded_targets.append(target_id)

        else:
            discarded_targets.append(target_id)

    return chosen_target_primers, discarded_targets
