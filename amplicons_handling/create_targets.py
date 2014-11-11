
import re
import hashlib
import csv
import argparse
from primer_design import primer3_design, bowtie2_design, sort_unique_primers
from primers_insertion import create_primers_in_db
from positioning import create_next_primers_plates, insertion_plates_to_db
from linapp.models import User, TargetEnrichmentType


if '__main__' == __name__:
    from django.core.management import setup_environ
    import sys
    sys.path.append('/home/ofirr/CLineage/')
    from clineage import settings
    setup_environ(settings)

from linapp.models import TargetType, Assembly, Target, Sequence, Chromosome, Microsatellite, SNP


columns_case_dict = {
                    ('Assembly', 'Chromosome', 'Start', 'End'):'Nameless',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Name'):'NoSec',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Name', 'Sequence'):'Plain',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Name', 'Sequence', 'SNP'):'SNP',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Name', 'Sequence', 'Repeat_Type', 'Repeat_Unit_Length', 'Repeat_length'):'MicroSatellite',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Partner'):'Nameless',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Partner', 'Name'):'NoSec',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Partner', 'Name', 'Sequence'):'Plain',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Partner', 'Name', 'Sequence', 'SNP'):'SNP',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Partner', 'Name', 'Sequence', 'Repeat_Type', 'Repeat_Unit_Length', 'Repeat_length'):'MicroSatellite'
                    }


case_target_type = {
                    'Nameless': TargetType.objects.get(name='NoSeq'),
                    'NoSec': TargetType.objects.get(name='NoSeq'),
                    'Plain': TargetType.objects.get(name='Plain'),
                    'SNP': TargetType.objects.get(name='SNP'),
                    'MicroSatellite': TargetType.objects.get(name='Microsatellite'),
}


def get_or_create_sequence(seq):
    if not re.match('^[ACTGactg]+$', seq.strip()):
        print 'unsupported characters in input sequence {}'.format(seq)
        raise
    try:
        sequence = Sequence.objects.get(hash=hashlib.md5(seq).hexdigest())
    except Sequence.DoesNotExist:
        sequence = Sequence(length=len(seq), sequence=seq, hash=hashlib.md5(seq).hexdigest())
        sequence.save()
    return sequence


def get_case_from_columns(columns):
    col_tup = tuple(columns)
    if columns_case_dict[col_tup]:
        return columns_case_dict[col_tup]
    print 'unsupported columns structure' #add helpful example
    raise


def clean_chromosome_name(raw_chr_name):
    if re.match('^[0-9]+$', raw_chr_name):
        return raw_chr_name
    if re.match('^[Cc]hr[0-9XYMxym]+$', raw_chr_name):
        return raw_chr_name[3:]
    if re.match('^[XYMxym]$', raw_chr_name):
        return raw_chr_name.upper()
    print 'unsupport chromosome name format'
    raise


def parse_commons(row_dict):
    assem = Assembly.objects.get(friendly_name=row_dict['Assembly'])
    chromosome_name = clean_chromosome_name(row_dict['Chromosome'])
    chrom = Chromosome.objects.get(name=chromosome_name, assembly=assem)
    start_pos = int(row_dict['Start'])
    end_pos = int(row_dict['End'])
    partner = None
    if row_dict['Partner']:
        try:
            partner = User.objects.get(username=row_dict['Partner'])
        except User.DoesNotExist:
            #warning
            pass
    return chrom, start_pos, end_pos, partner


def process_row(row_dict, case):
    chrom, start_pos, end_pos, partner = parse_commons(row_dict)
    name = '{}_{}_{}'.format(chrom.name, start_pos, end_pos)

    sequence = get_or_create_sequence(chrom.getdna(start_pos, end_pos))
    tgtype = case_target_type[case]

    if case in ['NoSec', 'Plain', 'SNP', 'MicroSatellite']:
        name = row_dict['Name']

    if case in ['Plain', 'SNP', 'MicroSatellite']:
        sequence = get_or_create_sequence(row_dict['Sequence'])
        start_pos, end_pos = chrom.locate(start_pos, end_pos, sequence.sequence)

    if case in ['SNP']:
        modified = row_dict['SNP']
        mutation = '{}>{}'.format(sequence.sequence, modified)
        obj, created = SNP.objects.get_or_create(
            start_pos=start_pos, end_pos=end_pos,
            defaults={'name': name,
                      'type': tgtype,
                      'chromosome': chrom,
                      'referencevalue': sequence,
                      'partner': partner,
                      'mutation':mutation,
                      'modified':modified}
        )
        return obj, created

    if case in ['MicroSatellite']:
        repeat_type = row_dict['Repeat_Type']
        repeat_unit_length = int(row_dict['Repeat_Unit_Length'])
        repeat_len = int(row_dict['Repeat_Length'])
        obj, created = Microsatellite.objects.get_or_create(
            start_pos=start_pos, end_pos=end_pos,
            defaults={'name':name,
                      'type':tgtype,
                      'chromosome':chr,
                      'referencevalue':sequence,
                      'partner': partner,
                      'repeat_unit_len':repeat_unit_length,
                      'repeat_unit_type':repeat_type,
                      'repeat_number':repeat_len}
        )
        return obj, created

    obj, created = Target.objects.get_or_create(
        start_pos=start_pos, end_pos=end_pos,
        defaults={'name':name,
                  'type':tgtype,
                  'chromosome':chr,
                  'referencevalue':sequence,
                  'partner': partner}
    )
    return obj, created


def proccess_input_target_file(input_file):
    obj_list = []
    with open(input_file, 'rb') as f:
        dialect = csv.Sniffer().sniff(f.read(1000))
        f.seek(0)
        rdr = csv.DictReader(f, dialect=dialect)
        row_case = get_case_from_columns(rdr.fieldnames)
        for row in rdr:
            obj, created = process_row(row, row_case)
            obj_list.appent(obj)
            print "Target {} INFO: {}".format(obj, created)
    return obj_list



if '__main__' == __name__:
    parser = argparse.ArgumentParser(description='Analyses hist-pairs file')
    parser.add_argument('-i', '--input', type=str, dest='input_file', help='path for target table file')
    parser.add_argument('-n', '--name', type=str, dest='input_name', help='path for target table file')
    parser.add_argument('-o', '--output', type=str, dest='output_name', help='output file name prefix for the ')
    parser.add_argument('-b', '--bowtie2Index', type=str, dest='bowtie2_index', help='path for target table file')
    parser.add_argument('-t', '--tails', type=bool, dest='tails', help='primers have tails?')
    args = parser.parse_args()
    input_file = args.input_file
    input_name = args.input_name
    output_name = args.output_name
    bowtie2_index = args.bowtie2_index
    xls_name = ("PrimerOrder{}.xls".format(str(output_name)))
    is_tails = args.tails
    no_tails_te_type, tails_te_type = TargetEnrichmentType.objects.all()
    if is_tails:
        te_type = tails_te_type
    else:
        te_type = no_tails_te_type
    obj_list = proccess_input_target_file(input_file)
    primer3_name_file = primer3_design(obj_list, input_name, output_name)
    sam_file, primer_data_check, target_primers = bowtie2_design(input_name, output_name, bowtie2_index, primer3_name_file)
    chosen_target_primers, discarded_targets = sort_unique_primers(sam_file, target_primers)
    create_primer_pairs = create_primers_in_db(chosen_target_primers, te_type)
    plate_united, plate_fw, plate_rev = create_next_primers_plates(obj_list[0].chromosome.assembly)
    insertion_plates_to_db(create_primer_pairs, xls_name, plate_united, plate_fw, plate_rev)



