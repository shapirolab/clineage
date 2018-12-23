from genomes.models import Chromosome, DNASlice
from targeted_enrichment.amplicons.models import PlainTargetedAmplicon, AmpliconCollection
from targeted_enrichment.planning.models import UGSPlus, UGSMinus, TargetEnrichment, Microsatellite
from targeted_enrichment.reagents.models import PCR1PrimerPairTER
from primers.synthesis.models import PCR1PlusPrimer, PCR1MinusPrimer

from utils.ms_utils import ms_type_permutations
from setup import irac1, irac2
import pandas as pd

def get_chromosome(row):
    chromosome_name = row['CHROM'] if 'chr' not in row['CHROM'] else row['CHROM'][3:]
    c = Chromosome.objects.get(name=chromosome_name)
    return c


def document_dna_slice(c, row):
    s, created = DNASlice.objects.get_or_create(
        chromosome=c,
        start_pos=int(row['START']),
        end_pos=int(row['END']),
    )
    return s


def document_microsatellite(s, row):
    repeat_unit_len = int(row['MOTIF_LEN'])
    repeat_unit_ref_seq = s.sequence[:repeat_unit_len]
    repeat_type_variants = set(ms_type_permutations(repeat_unit_ref_seq))
    canon = min(repeat_type_variants, key=lambda x: x.seq)
    ms, created = Microsatellite.objects.get_or_create(
        name=row['NAME'],
        slice=s,
        repeat_unit_len=row['MOTIF_LEN'],
        repeat_unit_type=canon.seq.decode(),
        repeat_number=row['NUM_COPIES'],
        repeat_unit_ref_seq=repeat_unit_ref_seq,
        planning_version=0,
    )
    return ms, created


def document_targeted_amplicon(c, row):
    amplicon_s, created = DNASlice.objects.get_or_create(
        chromosome=c,
        start_pos=int(row['START']) - 100,
        end_pos=int(row['END']) + 100,
    )
    left_s, created = DNASlice.objects.get_or_create(
        chromosome=c,
        start_pos=int(row['START']) - 100,
        end_pos=int(row['START']) - 80,
    )
    right_s, created = DNASlice.objects.get_or_create(
        chromosome=c,
        start_pos=int(row['END']) + 80,
        end_pos=int(row['END']) + 100,
    )
    left_ugs, created = UGSPlus.objects.get_or_create(slice=left_s)
    right_ugs, created = UGSMinus.objects.get_or_create(slice=right_s)
    left_primer, created = PCR1PlusPrimer.objects.get_or_create(
        name='',
        ugs=left_ugs,
        iraft=irac1,
    )
    right_primer, created = PCR1MinusPrimer.objects.get_or_create(
        name='',
        ugs=right_ugs,
        iraft=irac2,
    )
    te, created = TargetEnrichment.objects.get_or_create(
        chromosome=c,
        left=left_ugs,
        right=right_ugs,
        planning_version=1,
    )
    amplicon, created = PlainTargetedAmplicon.objects.get_or_create(
        slice=amplicon_s,
        left_ugs=left_ugs,
        right_ugs=right_ugs,
    )
    ter, created = PCR1PrimerPairTER.objects.get_or_create(
        te=te,
        passed_validation=False,
        left_primer=left_primer,
        right_primer=right_primer,
        amplicon=amplicon,
    )
    return amplicon, ter


def document_targets(regions_path):
    df = pd.read_csv(regions_path, sep='\t', names=["CHROM", "START", "END", "MOTIF_LEN", "NUM_COPIES", "NAME"])
    amplicons = []
    ters = []
    for i, row in df.iterrows():
        c = get_chromosome(row)
        s = document_dna_slice(c, row)
        ms, created = document_microsatellite(s, row)
        amplicon, ter = document_targeted_amplicon(c, row)
        amplicons.append(amplicon)
        ters.append(ter)
    return amplicons, ters
