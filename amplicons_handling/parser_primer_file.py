__author__ = 'veronika'
import re
import hashlib
from .primers_insertion import get_or_create_sequence, PrimerLocationError
from django.contrib.auth.models import User
from genomes.models import TargetType, Assembly, Chromosome, Microsatellite, SNP
from targeted_enrichment.planning.models import Target, Microsatellite, SNP
from utils.SequenceManipulations import complement

columns_case_dict = {
                    ('Assembly', 'Chromosome', 'Start', 'End'):'Nameless',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Name'):'NoSec',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Name', 'Sequence'):'Plain',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Name', 'Sequence', 'SNP'):'SNP',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Name', 'Sequence', 'Repeat_Type', 'Repeat_Unit_Length', 'Repeat_Length'):'MicroSatellite',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Partner'):'Nameless',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Partner', 'Name'):'NoSec',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Partner', 'Name', 'Sequence'):'Plain',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Partner', 'Name', 'Sequence', 'SNP'):'SNP',
                    ('Assembly', 'Chromosome', 'Start', 'End', 'Partner', 'Name', 'Sequence', 'Repeat_Type', 'Repeat_Unit_Length', 'Repeat_Length'):'MicroSatellite'
                    }


case_target_type = {
                    'Nameless': TargetType.objects.get(name='NoSeq'),
                    'NoSec': TargetType.objects.get(name='NoSeq'),
                    'Plain': TargetType.objects.get(name='Plain'),
                    'SNP': TargetType.objects.get(name='SNP'),
                    'MicroSatellite': TargetType.objects.get(name='Microsatellite'),
}


def get_case_from_columns(columns):
    col_tup = tuple(columns)
    if columns_case_dict[col_tup]:
        return columns_case_dict[col_tup]
    print('unsupported columns structure') #add helpful example
    raise


def clean_chromosome_name(raw_chr_name):
    if re.match('^[0-9]+$', raw_chr_name):
        return raw_chr_name
    if re.match('^[Cc]hr[0-9XYMxym]+$', raw_chr_name):
        return raw_chr_name[3:]
    if re.match('^[XYMxym]$', raw_chr_name):
        return raw_chr_name.upper()
    print('unsupport chromosome name format')
    raise


def parse_commons(row_dict):
    assem = Assembly.objects.get(friendly_name=row_dict['Assembly'])
    chromosome_name = clean_chromosome_name(row_dict['Chromosome'])
    chrom = Chromosome.objects.get(name=chromosome_name, assembly=assem)
    start_pos = int(row_dict['Start'])
    end_pos = int(row_dict['End'])
    partner = None
    if 'Partner' in row_dict:
        try:
            partner = User.objects.get(username=row_dict['Partner'])
        except User.DoesNotExist:
            #warning
            partner = None
    return chrom, start_pos, end_pos, partner


def snp_object(row_dict, sequence, start_pos, end_pos, name, tgtype, chrom, partner):
    modified = row_dict['SNP']
    mutation = '{}>{}'.format(sequence.sequence, modified)
    obj, created = SNP.objects.get_or_create(
        start_pos=start_pos, end_pos=end_pos,
        defaults={'name': name,
                  'type': tgtype,
                  'chromosome': chrom,
                  'referencevalue': sequence,
                  'mutation':mutation,
                  'modified':modified}
    )
    if not partner or partner in obj.partner.all():
        return obj, created
    obj.partner.add(partner)
    obj.save()
    return obj, created


def locate_sequence_on_strand(chrom, start_pos, end_pos, sequence, margins):
    try:
        ms_s, ms_e = chrom.locate(start_pos,
                                  end_pos,
                                  sequence.sequence,
                                  padding=margins)
        if ms_s != start_pos or ms_e != end_pos:
            print('WARN: input indexes are off and were corrected')
    except ValueError:
        try:
            ms_s, ms_e = chrom.locate(start_pos,
                                      end_pos,
                                      complement(sequence.sequence)[::-1],
                                      padding=margins)
            if ms_s != start_pos or ms_e != end_pos:
                print('WARN: input indexes are off and were corrected')
        except ValueError:
            print(chrom.name, chrom.assembly, start_pos, end_pos, sequence.sequence)
            raise PrimerLocationError
    return ms_s, ms_e


def microsatellite_object(row_dict, sequence, start_pos, end_pos, name, tgtype, chrom, partner):
    repeat_type = row_dict['Repeat_Type']
    repeat_unit_length = int(row_dict['Repeat_Unit_Length'])
    repeat_len = int(row_dict['Repeat_Length'])
    ######################################
    overlapping = False
    if Microsatellite.objects.filter(chromosome=chrom)\
                             .filter(start_pos__lte=start_pos)\
                             .filter(end_pos__gte=start_pos) or \
        Microsatellite.objects.filter(chromosome=chrom)\
                              .filter(start_pos__lte=end_pos)\
                              .filter(end_pos__gte=end_pos):
        overlapping = True
    ######################################
    obj, created = Microsatellite.objects.get_or_create(
        start_pos=start_pos, end_pos=end_pos,
        defaults={'name': name,
                  'type': tgtype,
                  'chromosome': chrom,
                  'referencevalue': sequence,
                  'repeat_unit_len': repeat_unit_length,
                  'repeat_unit_type': repeat_type,
                  'repeat_number': repeat_len}
    )
    if overlapping and created:
        print('WARN: overlapping MS for ', start_pos, end_pos, name, tgtype, chrom)
    if not partner or partner in obj.partner.all():
        return obj, created
    obj.partner.add(partner)
    obj.save()
    return obj, created


def nosec_object(sequence, start_pos, end_pos, name, tgtype, chrom, partner):
    obj, created = Target.objects.get_or_create(
        start_pos=start_pos, end_pos=end_pos, chromosome=chrom,
        defaults={'name': name,
                  'type': tgtype,
                  'referencevalue': sequence}
        )
    if not partner or partner in obj.partner.all():
        return obj, created
    obj.partner.add(partner)
    obj.save()
    return obj, created


def process_row(row_dict, case, margins=10):
    chrom, start_pos, end_pos, partner = parse_commons(row_dict)
    name = '{}_{}_{}'.format(chrom.name, start_pos, end_pos)

    sequence = get_or_create_sequence(chrom.getdna(start_pos, end_pos))
    tgtype = case_target_type[case]

    if case in ['NoSec', 'Plain', 'SNP', 'MicroSatellite']:
        name = row_dict['Name']

    if case in ['Plain', 'SNP', 'MicroSatellite']:
        sequence = get_or_create_sequence(row_dict['Sequence'])
        start_pos, end_pos = locate_sequence_on_strand(chrom, start_pos, end_pos, sequence, margins)

    if case in ['SNP']:
        return snp_object(row_dict, sequence, start_pos, end_pos, name, tgtype, chrom, partner)

    if case in ['MicroSatellite']:
        return microsatellite_object(row_dict, sequence, start_pos, end_pos, name, tgtype, chrom, partner)

    return nosec_object(sequence, start_pos, end_pos, name, tgtype, chrom, partner)
