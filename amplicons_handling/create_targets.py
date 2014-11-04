
import re
import hashlib
import csv
import argparse

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
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Name', 'Sequence', 'Repeat_Type', 'Repeat_Unit_Length', 'Repeat_length'):'MicroSatellite'
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
    if columns_case_dict[columns]:
        return columns_case_dict[columns]
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
    chr = Chromosome.objects.get(name=chromosome_name, assembly=assem)
    start_pos = int(row_dict['Start'])
    end_pos = int(row_dict['End'])
    return chr, start_pos, end_pos


def process_row(row_dict, case):
    chrom, start_pos, end_pos = parse_commons(row_dict)
    name = '{}_{}_{}'.format(chrom.name, start_pos, end_pos)
    sequence = get_or_create_sequence(chrom.getdna(start_pos, end_pos))
    tgtype = case_target_type[case]
    if case in ['NoSec', 'Plain', 'SNP', 'MicroSatellite']:
        name = row_dict['Name']

    if case in ['Plain', 'SNP', 'MicroSatellite']:
        sequence = get_or_create_sequence(row_dict['Sequence'])
        start_pos, end_pos = chrom.locate(start_pos, end_pos, sequence)

    if case in ['SNP']:
        modified = row_dict['SNP']
        mutation = '{}>{}'.format(sequence, modified)
        obj, created = SNP.objects.get_or_create(start_pos=start_pos, end_pos=end_pos,
                                            defaults={'name': name, 'type': tgtype, 'chromosome': chrom, 'referencevalue': sequence, 'mutation':mutation, 'modified':modified})
        return obj, created

    if case in ['MicroSatellite']:
        repeat_type = row_dict['Repeat_Type']
        repeat_unit_length = int(row_dict['Repeat_Unit_Length'])
        repeat_len = int(row_dict['Repeat_Length'])
        obj, created = Microsatellite.objects.get_or_create(start_pos=start_pos, end_pos=end_pos,
                                            defaults={'name':name, 'type':tgtype, 'chromosome':chr, 'referencevalue':sequence, 'repeat_unit_len':repeat_unit_length, 'repeat_unit_type':repeat_type, 'repeat_number':repeat_len})
        return obj, created

    obj, created = Target.objects.get_or_create(start_pos=start_pos, end_pos=end_pos,
                                            defaults={'name':name, 'type':tgtype, 'chromosome':chr, 'referencevalue':sequence})
    return obj, created

if '__main__' == __name__:
    parser = argparse.ArgumentParser(description='Analyses hist-pairs file')
    parser.add_argument('-i', '--input', type=str, dest='input_file', help='path for target table file')
    args = parser.parse_args()
    input_file = args.input_file
    with open(input_file, 'rb') as f:
        dialect = csv.Sniffer().sniff(f.read(1000))
        f.seek(0)
        rdr = csv.reader(f, dialect=dialect)
        case = get_case_from_columns(rdr.fieldnames)
        for row in rdr:
            process_row(row, case)