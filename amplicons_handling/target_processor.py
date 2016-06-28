import csv
import argparse
from .parser_primer_file import get_case_from_columns, process_row
from .primer_design import primer3_design, bowtie2_design, sort_best_primers
from .primers_insertion import create_primers_in_db
from .positioning import insertion_plates_to_db, create_primer_order_file_xls
from django.core.management import setup_environ
import sys
sys.path.append('/home/ofirr/CLineage')
from clineage import settings
from genomes.models import TargetEnrichmentType, PrimerTail
# setup_environ(settings)


def create_target_in_db(obj_list,
                        input_name,
                        output_name,
                        min_size,
                        max_size,
                        primer_num_rerun,
                        margins):
    primer3_name_file = primer3_design(obj_list,
                                       input_name,
                                       output_name,
                                       min_size,
                                       max_size,
                                       primer_num_rerun,
                                       margins)
    sam_file, primer_data_check, target_primers = bowtie2_design(input_name, output_name, bowtie2_index, primer3_name_file)
    chosen_target_primers, discarded_targets = sort_best_primers(sam_file, target_primers, margins=max_size)

    ptf, ptr = PrimerTail.objects.all()
    te_list, colliding_amplicons = create_primers_in_db(chosen_target_primers, te_type, pf_tail=ptf, pr_tail=ptr)
    pairs_plates, stk_fw_plates, stk_rv_plates = insertion_plates_to_db(te_list, assembly)
    create_primer_order_file_xls(stk_fw_plates, stk_rv_plates, xls_name)