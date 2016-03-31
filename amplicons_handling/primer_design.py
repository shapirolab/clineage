__author__ = 'veronika'

import os
import re
import argparse
import pysam
from django.conf import settings
from collections import defaultdict
from collections import Counter
from frogress import bar
from .primers_insertion import create_primers_in_db
from .positioning import insertion_plates_to_db, create_primer_order_file_xls
from .primers_insertion import check_primers, AmpliconCollisionError, PrimerLocationError
from genomes.models import TargetEnrichmentType, PrimerTail
from targeted_enrichment.planning.models import Target


def create_amplicons_for_primer3(target):
    amplicon = target.chromosome.getdna(target.start_pos-1000, target.end_pos+1000)

    return amplicon


def primer3_design(obj_list, input_name, output_name,
                   min_size=130, max_size=400, primer_num_rerun=10000, margins=0):
    primer3_input = ("{}.txt".format(str(input_name)))
    with open(primer3_input, 'w+') as primer3_file:
        primer3_file.write('PRIMER_TASK=pick_detection_primers\nPRIMER_OPT_SIZE=23\nPRIMER_MIN_SIZE=20\n'
                                    'PRIMER_MAX_SIZE=27\nPRIMER_PRODUCT_SIZE_RANGE={}-{}\nP3_FILE_FLAG=0\n'
                                    'PRIMER_EXPLAIN_FLAG=1\nPRIMER_MIN_TM=51\nPRIMER_OPT_TM=55\nPRIMER_MAX_TM=60\n'
                                    'PRIMER_SALT_CORRECTIONS=1\nPRIMER_TM_FORMULA=1\nPRIMER_PAIR_MAX_DIFF_TM=3\n'
                                    'PRIMER_NUM_RETURN={}\nPRIMER_FILE_FLAG=0\n'.format(min_size,
                                                                                        max_size,
                                                                                        primer_num_rerun))
        for target in obj_list:
            amplicon = create_amplicons_for_primer3(target)
            primer3_file.write('SEQUENCE_ID={}\nSEQUENCE_TEMPLATE={}\n'
                               'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST={},{},{},{}\n=\n'.format(target.id,
                                                                                             amplicon,
                                                                                             1000-(max_size-(len(target.referencevalue.sequence)/2+margins)),
                                                                                             max_size-(len(target.referencevalue.sequence)/2+margins),
                                                                                             1000+len(target.referencevalue.sequence),
                                                                                             max_size-(len(target.referencevalue.sequence)/2+margins)))

    # Run the primer3 on the input
    primer3_output = ("{}_primer3.txt".format(str(output_name)))
    s = '{} < {} > {}'.format(settings.PRIMER3_PATH, primer3_input, primer3_output)
    os.system(s)
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
            target_primers[seq_id]['LEFT'][i] = primer_left
            target_primers[seq_id]['RIGHT'][i] = primer_right
    return target_primers


def bowtie2_design(input_fasta_file, output_file, bowtie2_index, output_name):
    target_primers = parse_primer3_output(output_name)
    out_string = ''
    for seq_id in list(target_primers.keys()):
        for primer_number in list(target_primers[seq_id]['LEFT'].keys()):
            out_string += '>PRIMER_LEFT_{}_{}\n{}\n'.format(primer_number, seq_id, target_primers[seq_id]['LEFT'][primer_number])
            out_string += '>PRIMER_RIGHT_{}_{}\n{}\n'.format(primer_number, seq_id, target_primers[seq_id]['RIGHT'][primer_number])
    primer_data_check = '{}.fa'.format(str(input_fasta_file))
    print('writing primers fasta file {}'.format(primer_data_check))
    with open(primer_data_check, 'w+') as primers_output:
        primers_output.write(out_string)
    sam_file = '{}.sam'.format(str(output_file))
    print('running bowtie sam file output:{}'.format(sam_file))
    s = '{} -k 2 {} -f -U {} -S {}'.format(settings.BOWTIE2_PATH, bowtie2_index, primer_data_check, sam_file)
    os.system(s)
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


def sort_unique_primers(sam_file, target_primers, margins=300):
    name_count, primers_names = primer_count_from_sam_file(sam_file)
    chosen_target_primers = defaultdict(lambda: defaultdict(str))
    discarded_targets = []
    for target_id in bar(list(target_primers.keys())):
        left_primer = None
        for primer_number_fr in sorted(target_primers[target_id]['LEFT'].keys()):
            bowtie2_key_left = 'PRIMER_LEFT_{}_{}'.format(primer_number_fr, target_id)
            if name_count[bowtie2_key_left] == 1:
                left_primer = target_primers[target_id]['LEFT'][primer_number_fr]
                right_primer = None
                for primer_number_rv in sorted(target_primers[target_id]['RIGHT'].keys()):
                    bowtie2_key_right = 'PRIMER_RIGHT_{}_{}'.format(primer_number_rv, target_id)
                    if name_count[bowtie2_key_right] == 1:
                        right_primer = target_primers[target_id]['RIGHT'][primer_number_rv]
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


# if '__main__' == __name__:
#     parser = argparse.ArgumentParser(description='Analyses hist-pairs file')
#     parser.add_argument('-i', '--input', type=str, dest='obj_list', help='list af targets for primers design')
#     parser.add_argument('-p', '--inputFolder', type=str, dest='input_folder', help='path for target folder')
#     parser.add_argument('-n', '--name', type=str, dest='input_name', help='files prefix for the targets')
#     parser.add_argument('-o', '--output', type=str, dest='output_name', help='output file name prefix for the targets')
#     parser.add_argument('-b', '--bowtie2Index', type=str, dest='bowtie2_index', help='path for bowtie2_index files')
#     parser.add_argument('-t', '--tails', type=bool, dest='tails', help='primers have tails?')
#     parser.add_argument('-r', '--primerNumReRun', type=int, dest='primer_num_rerun', default=10000, help='number of pair primers to constract')
#     parser.add_argument('-m', '--margins', type=int, dest='margins', default=0, help='margins for the desired sequencing region')
#     parser.add_argument('-a', '--minSize', type=int, dest='min_size', default=130, help='minimus amplicon size')
#     parser.add_argument('-z', '--maxSize', type=int, dest='max_size', default=300, help='maximum amplicon size')
#     parser.add_argument('-s', '--inSilico', type=bool, dest='insilico_test', default=True, help='Run in silico test for primers')
#     args = parser.parse_args()
#     obj_list = args.obj_list
#     input_file = args.input_file
#     input_folder = args.input_folder
#     input_name = args.input_name
#     output_name = args.output_name
#     bowtie2_index = args.bowtie2_index
#     primer_num_rerun = args.primer_num_rerun
#     margins = args.margins
#     min_size = args.min_size
#     max_size = args.max_size
#     insilico_test = args.insilico_test
#     xls_name = ("PrimerOrder{}.xls".format(str(output_name)))
#     is_tails = args.tails
#     no_tails_te_type, tails_te_type = TargetEnrichmentType.objects.all()
#     if is_tails:
#         te_type = tails_te_type
#     else:
#         te_type = no_tails_te_type
#     primer3_name_file = primer3_design(obj_list,
#                                        input_name,
#                                        output_name,
#                                        min_size,
#                                        max_size,
#                                        primer_num_rerun,
#                                        margins)
#     sam_file, primer_data_check, target_primers = bowtie2_design(input_name, output_name, bowtie2_index, primer3_name_file)
#     chosen_target_primers, discarded_targets = sort_unique_primers(sam_file, target_primers, margins=max_size)
#     print 'amount of chosen tragets: {}, amount of discarded targets: {}'.format(len(chosen_target_primers), len(discarded_targets))
#     ptf, ptr = PrimerTail.objects.all()
#     te_list = create_primers_in_db(chosen_target_primers, te_type, margins, insilico_test, pf_tail=ptf, pr_tail=ptr)
#     pairs_plates, stk_fw_plates, stk_rv_plates = insertion_plates_to_db(te_list)
#     create_primer_order_file_xls(stk_fw_plates, stk_rv_plates, xls_name)
