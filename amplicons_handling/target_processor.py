import csv
import argparse
from .primer_design import primer3_design, sort_best_primers
from .primers_insertion import create_target_enrichment_in_db
from .positioning import insertion_OM_to_db, create_primer_order_file_xls
from primers.parts.models import IlluminaReadingAdaptor1ForTail, IlluminaReadingAdaptor2ForTail


def create_target_in_db(target_list, om_name, xls_name,  **kwargs):

    tate_tuple = []
    ira1ft = IlluminaReadingAdaptor1ForTail.objects.get(ira__name='Illumina Standard Reading Adaptor1', tail_length=27)
    ira2ft = IlluminaReadingAdaptor2ForTail.objects.get(ira__name='Illumina Standard Reading Adaptor2', tail_length=27)

    for target in target_list:
        primer3_output = primer3_design(target, **kwargs)
        chosen_target_primers, discarded_targets = sort_best_primers(primer3_output, **kwargs)

        if chosen_target_primers:
            # create TE
            ta, te = create_target_enrichment_in_db(chosen_target_primers)
            tate_tuple.append((ta, te))
    print('amount of TE created: {}, amount of all targets submitted: {}'.format(len(tate_tuple), len(target_list)))
    om_mix = insertion_OM_to_db(tate_tuple, om_name, ira1ft, ira2ft)
    create_primer_order_file_xls(om_mix, xls_name)
