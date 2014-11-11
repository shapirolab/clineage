__author__ = 'veronika'
from linapp.models import Target
import os
import re
import pysam
from collections import defaultdict
from collections import Counter


# if '__main__' == __name__:
#     from django.core.management import setup_environ
#     import sys
#     sys.path.append('/home/ofirr/CLineage/')
#     from clineage import settings
#     setup_environ(settings)


def create_amplicons_for_primer3(target, margins):
    amplicon = target.chromosome.getdna(target.start_pos-margins, target.end_pos+margins)

    return amplicon


def primer3_design(obj_list, input_name, output_name, primer_num_rerun=10000, margins=200, seq_region=100):
    primer3_input = ("{}.txt".format(str(input_name)))
    with open(primer3_input, 'w+') as primer3_file:
        primer3_file.write('PRIMER_TASK=pick_detection_primers\nPRIMER_OPT_SIZE=23\nPRIMER_MIN_SIZE=20\n'
                                    'PRIMER_MAX_SIZE=27\nPRIMER_PRODUCT_SIZE_RANGE=30-1000\nP3_FILE_FLAG=0\n'
                                    'PRIMER_EXPLAIN_FLAG=1\nPRIMER_MIN_TM=51\nPRIMER_OPT_TM=55\nPRIMER_MAX_TM=60\n'
                                    'PRIMER_SALT_CORRECTIONS=1\nPRIMER_TM_FORMULA=1\nPRIMER_PAIR_MAX_DIFF_TM=3\n'
                                    'PRIMER_NUM_RETURN={}\nPRIMER_FILE_FLAG=0\n').format(primer_num_rerun)
        for target in obj_list:
            amplicon = create_amplicons_for_primer3(target, margins)
            primer3_file.write('SEQUENCE_ID={}\nSEQUENCE_TEMPLATE={}\n'
                               'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,{},{},{}\n=\n'.format(target.id,
                                                                                            amplicon,
                                                                                            seq_region,
                                                                                            len(amplicon)-seq_region,
                                                                                            seq_region))

    #Run the primer3 on the input
    primer3_output = ("{}.txt".format(str(output_name)))
    s = '/net/mraid11/export/data/dcsoft/home/Adam/PCRPrimersDesign/Primer3/primer3_core < {} > {}'.format(primer3_input,primer3_output)
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
    for seq_id in target_primers.keys():
        for primer_number in enumerate(target_primers[seq_id]['LEFT'].keys()):
            out_string += '>PRIMER_LEFT_{}_{}\n{}\n'.format(primer_number, seq_id, target_primers[seq_id]['LEFT'][primer_number])
            out_string += '>PRIMER_RIGHT_{}_{}\n{}\n'.format(primer_number, seq_id, target_primers[seq_id]['RIGHT'][primer_number])
        primer_data_check = '{}.fa'.format(str(input_fasta_file))
        primers_output = open(primer_data_check, 'w+')
        primers_output.write(out_string)
        primers_output.close()
        sam_file = '{}.sam'.format(str(output_file))
        s = '/net/mraid11/export/data/dcsoft/home/Adam/Software/bowtie2-2.2.2/bowtie2 -k 2 {} -f -U {} -S {}'.format(bowtie2_index, primer_data_check, sam_file)
        os.system(s)
        return sam_file, primer_data_check, target_primers


def primer_count_from_sam_file(sam_file):
    samfile = pysam.Samfile(sam_file, "r")
    primers_names = []
    for sam_primer in samfile:
        primers_names.append(sam_primer.qname)
    name_count = Counter(primers_names)
    return name_count


def sort_unique_primers(sam_file, target_primers):
    name_count = primer_count_from_sam_file(sam_file)
    chosen_target_primers = defaultdict(lambda: defaultdict(str))
    discarded_targets = []
    for target_id in target_primers.keys():
        #left
        left_primer = None
        for primer_number in sorted(target_primers[target_id]['LEFT'].keys()):
            bowtie2_key = 'PRIMER_LEFT_{}_{}'.format(primer_number, target_id)
            if name_count[bowtie2_key] == 1:
                left_primer = target_primers[target_id]['LEFT'][primer_number]
                break
        #right
        right_primer = None
        for primer_number in sorted(target_primers[target_id]['LEFT'].keys()):
            bowtie2_key = 'PRIMER_RIGHT_{}_{}'.format(primer_number, target_id)
            if name_count[bowtie2_key] == 1:
                right_primer = target_primers[target_id]['RIGHT'][primer_number]
                break

        if left_primer and right_primer:
            chosen_target_primers[target_id]['LEFT'] = left_primer
            chosen_target_primers[target_id]['RIGHT'] = right_primer
        else:
            discarded_targets.append(target_id)

    return chosen_target_primers, discarded_targets



# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser.add_argument("loci_names_file", help=".bed file name for generating primers",
#                         type=str)
#     parser.add_argument("ListFileName", help="prefix name for all the files that are generated",
#                         type=str)
#     args = parser.parse_args()
#
#     create_primers_for_new_targets(args.loci_names_file, args.ListFileName)