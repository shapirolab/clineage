import csv
import argparse
from parser_primer_file import get_case_from_columns, process_row
from primer_design import primer3_design, bowtie2_design, sort_unique_primers
from primers_insertion import create_primers_in_db
from positioning import insertion_plates_to_db, create_primer_order_file_xls
from django.core.management import setup_environ
import sys
sys.path.append('/home/ofirr/CLineage')
from clineage import settings
from genomes.models import TargetEnrichmentType, PrimerTail
setup_environ(settings)


def proccess_input_target_file(input_file, margins):
    obj_list = []
    assem = None
    with open(input_file, 'rb') as f:
        dialect = csv.Sniffer().sniff(f.read(10000))
        f.seek(0)
        rdr = csv.DictReader(f, dialect=dialect)
        row_case = get_case_from_columns(rdr.fieldnames)
        for row in rdr:
            if not assem:
                assem = row['Assembly']
            obj, created = process_row(row, row_case, margins)
            obj_list.append(obj)
            print "Target {} INFO: {}".format(obj, created)
    return obj_list, assem


if '__main__' == __name__:
    parser = argparse.ArgumentParser(description='Analyses hist-pairs file')
    parser.add_argument('-i', '--input', type=str, dest='input_file', help='path for target table file')
    parser.add_argument('-p', '--inputFolder', type=str, dest='input_folder', help='path for target folder')
    parser.add_argument('-n', '--name', type=str, dest='input_name', help='files prefix for the targets')
    parser.add_argument('-o', '--output', type=str, dest='output_name', help='output file name prefix for the targets')
    parser.add_argument('-b', '--bowtie2Index', type=str, dest='bowtie2_index', help='path for bowtie2_index files')
    parser.add_argument('-t', '--tails', type=bool, dest='tails', help='primers have tails?')
    parser.add_argument('-r', '--primerNumReRun', type=int, dest='primer_num_rerun', default=10000, help='number of pair primers to constract')
    parser.add_argument('-m', '--margins', type=int, dest='margins', default=0, help='margins for the desired sequencing region')
    parser.add_argument('-a', '--minSize', type=int, dest='min_size', default=130, help='minimus amplicon size')
    parser.add_argument('-z', '--maxSize', type=int, dest='max_size', default=400, help='maximum amplicon size')
    parser.add_argument('-s', '--inSilico', type=bool, dest='insilico_test', default=True, help='Run in silico test for primers')
    args = parser.parse_args()
    input_file = args.input_file
    input_folder = args.input_folder
    input_name = args.input_name
    output_name = args.output_name
    bowtie2_index = args.bowtie2_index
    primer_num_rerun = args.primer_num_rerun
    margins = args.margins
    min_size = args.min_size
    max_size = args.max_size
    insilico_test = args.insilico_test
    xls_name = ("PrimerOrder{}.xls".format(str(output_name)))
    is_tails = args.tails
    no_tails_te_type = TargetEnrichmentType.objects.get(name='PCR')
    tails_te_type = TargetEnrichmentType.objects.get(name='PCR_with_tails')
    if is_tails:
        te_type = tails_te_type
    else:
        te_type = no_tails_te_type
    obj_list, assembly = proccess_input_target_file(input_file, margins)
    primer3_name_file = primer3_design(obj_list,
                                       input_name,
                                       output_name,
                                       min_size,
                                       max_size,
                                       primer_num_rerun,
                                       margins)
    sam_file, primer_data_check, target_primers = bowtie2_design(input_name, output_name, bowtie2_index, primer3_name_file)
    chosen_target_primers, discarded_targets = sort_unique_primers(sam_file, target_primers, margins=max_size)
    print 'amount of chosen tragets: {}, amount of discarded targets: {}'.format(len(chosen_target_primers), len(discarded_targets))
    ptf, ptr = PrimerTail.objects.all()
    te_list, colliding_amplicons = create_primers_in_db(chosen_target_primers, te_type, pf_tail=ptf, pr_tail=ptr)
    pairs_plates, stk_fw_plates, stk_rv_plates = insertion_plates_to_db(te_list, assembly)
    create_primer_order_file_xls(stk_fw_plates, stk_rv_plates, xls_name)



