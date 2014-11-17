__author__ = 'veronika'

import os
import subprocess
import shutil
import csv
import re
from os import path
from collections import defaultdict
# from subprocess import getstatusoutput


# def run_command(cmd):
#     #print 'Running: "%s"' % cmd
#     status, text = subprocess.call(cmd)
#     exit_code = status >> 8 # high byte
#     signal_num = status % 256 # low byte
#     return text
import argparse

def build_primer3_input(loci_names_file,ListFileName):
    primer3_input = ("Primer3_Input_{}.txt".format(str(ListFileName)))
    #Prepare the input file for primer3
    if '.bed' in loci_names_file:
        short_name = loci_names_file.split(".")[0]
        fasta_primers = '{}.fasta'.format(str(short_name))
        # fasta_primers_file = open(fasta_primers, 'w+')
        system_call = ()

        s = 'bedtools getfasta -fi {} -bed {} -fo {}'.format(chrom.get_abs_path(), loci_names_file, fasta_primers)
        # s = 'bedtools getfasta -fi ~/Adam/GenomesData/Human/hg19/chrM.fa -bed {} -fo {}'.format(loci_names_file, fasta_primers)
        os.system(s)

    with open(fasta_primers, 'rb') as fasta_primers_file:
        if fasta_primers_file:
            if 'chr' in fasta_primers_file.readline():
                primer_names = []
                with open(loci_names_file,'rb') as loci_names:
                    loci_name_ind = csv.reader(loci_names, dialect='excel-tab')
                    for row in loci_name_ind:
                        primer_names.append(row[3])
                    #print loci_names.split('\t')[3]
                with open(primer3_input, 'w+') as primer3_file:
                    primer3_file.write('PRIMER_TASK=pick_detection_primers\nPRIMER_OPT_SIZE=23\nPRIMER_MIN_SIZE=20\n'
                                   'PRIMER_MAX_SIZE=27\nPRIMER_PRODUCT_SIZE_RANGE=30-1000\nP3_FILE_FLAG=0\nPRIMER_EXPLAIN_FLAG=1'
                                   '\nPRIMER_MIN_TM=51\nPRIMER_OPT_TM=55\nPRIMER_MAX_TM=60\nPRIMER_SALT_CORRECTIONS=1\n'
                                   'PRIMER_TM_FORMULA=1\nPRIMER_PAIR_MAX_DIFF_TM=3\nPRIMER_NUM_RETURN=10000\nPRIMER_FILE_FLAG=0\n')

                    i=0
                    for i, line in enumerate(fasta_primers_file):
                        if i%2 == 1:
                            continue
                        # name_for_loci = primer_names[i].split("\t")[3].strip()
                        target_for_loci = line
                        primer3_file.write('SEQUENCE_ID={}\nSEQUENCE_TEMPLATE={}SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,100,200,100\n=\n'.format(primer_names[i/2],target_for_loci))
                        i=i+1

    #Run the primer3 on the input
    primer3_output = ("Primer3_Output_{}.txt".format(str(ListFileName)))
    # primers_output_file = open(primer3_output, 'w+')
    s= '/net/mraid11/export/data/dcsoft/home/Adam/PCRPrimersDesign/Primer3/primer3_core < {} > {}'.format(primer3_input,primer3_output)
    os.system(s)

    #Go over the results of primer3 and parse all the primers candidates and Generate fasta file for primers alignment
    primer_data_check = 'Primers_availability_{}.fa'.format(str(ListFileName))

    fasta_output_file = open(primer3_output, 'rb')
    s = fasta_output_file.read()
    out_string = ''
    fasta_string = s.split('\n=\n')
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
            out_string += '>PRIMER_LEFT_{}_{}\n{}\n'.format(i, seq_id, primer_left)
            out_string += '>PRIMER_RIGHT_{}_{}\n{}\n'.format(i, seq_id, primer_right)

    fasta_output_file.close()
    primers_output = open(primer_data_check, 'w+')
    primers_output.write(out_string)
    primers_output.close()

    sam_file = '{}.sam'.format(str(ListFileName))
    s= '/net/mraid11/export/data/dcsoft/home/Adam/Software/bowtie2-2.2.2/bowtie2 -k 2 ~/Adam/GenomesData/Human/hg19/bowtie2_index/hg19 -f -U {} -S {}'.format(primer_data_check,sam_file)
    os.system(s)
    uniq_sorted = '{}_uniq_sorted.txt'.format(str(ListFileName))
    s="""cat {} | grep chr | awk '{{print $1}}' | sort | uniq -c | sort -n > {}""".format(sam_file,uniq_sorted)
    os.system(s)
    all_candidates = '{}_all_candidates.txt'.format(str(ListFileName))
    s="""cat {} | grep '>' | sed s/">"// > {}""".format(primer_data_check,all_candidates)
    os.system(s)
    primers_remove = '{}_primers_remove.txt'.format(str(ListFileName))
    s="""cat {} | awk '{{if ($1>=2) print $2}}' > {}""".format(uniq_sorted,primers_remove)
    os.system(s)
    with open(all_candidates,'rb') as all_candidates_file:
        all_candidates_names = set(all_candidates_file.read().split('\n'))
    with open(primers_remove,'rb') as primers_remove_file:
        primers_remove_names = set(primers_remove_file.read().split('\n'))
    valid_choices = all_candidates_names - primers_remove_names
    valid_choices = [name.split('_')[1:] for name in valid_choices] #row = [direction, ord, seqid]
    choices = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for row in valid_choices:
        choices[row[2]][row[0]][int(row[1])] = 1
    to_order = defaultdict(lambda: defaultdict(str))
    bad_primer_names = set(primer_names)
    order_primers_file = '{}_order_primers.txt'.format(str(ListFileName))
    order_primers_output = open(order_primers_file, 'w+')
    for seq_id in target_primers.keys():
        found = 0
        for i in range(10000):
            if not found:
                if choices[seq_id]['LEFT'][i]:
                    for j in range(10000):
                        if choices[seq_id]['RIGHT'][j]:
                            order_primers_output.write('{}\t{}\t{}\n'.format(str(seq_id), 'LEFT', str(target_primers[seq_id]['LEFT'][i])))
                            order_primers_output.write('{}\t{}\t{}\n'.format(str(seq_id), 'RIGHT', str(target_primers[seq_id]['RIGHT'][j])))
                            bad_primer_names = bad_primer_names - set([seq_id])
                            found = 1
                            break
    order_primers_output.close()
    bad_primers_names = '{}_bad_primers_names.txt'.format(str(ListFileName))
    primers_output = open(bad_primers_names, 'w+')
    primers_output.write('\n'.join(list(bad_primer_names)))
    primers_output.close()
# Now, run the alignment on all the candidates and sort them by running this code:



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("loci_names_file", help=".bed file name for generating primers",
                        type=str)
    parser.add_argument("ListFileName", help="prefix name for all the files that are generated",
                        type=str)
    args = parser.parse_args()

    build_primer3_input(args.loci_names_file, args.ListFileName)