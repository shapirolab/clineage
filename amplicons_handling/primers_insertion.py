__author__ = 'veronika'

from targeted_enrichment.planning.models import Target, TargetEnrichment, UGSPlus, UGSMinus
from targeted_enrichment.amplicons.models import UMITargetedAmplicon
from genomes.models import DNASlice


def create_ugs_fwd(chr, start_pos, end_pos):
    slice, c = DNASlice.objects.get_or_create(
        chromosome=chr,
        start_pos=start_pos,
        end_pos=end_pos,
    )
    slice.cache()

    ugs, c = UGSPlus.objects.get_or_create(
        slice=slice,
    )

    return ugs


def create_ugs_rev(chr, start_pos, end_pos):
    slice, c = DNASlice.objects.get_or_create(
        chromosome=chr,
        start_pos=start_pos,
        end_pos=end_pos,
    )
    slice.cache()

    ugs, c = UGSMinus.objects.get_or_create(
        slice=slice,
    )

    return ugs


def create_new_clice(chr, start_pos, end_pos):
    slice, c = DNASlice.objects.get_or_create(
        chromosome=chr,
        start_pos=start_pos,
        end_pos=end_pos,
    )
    slice.cache()
    return slice


def create_target_enrichment_in_db(chosen_target_primers):

    target = Target.objects.get(id=chosen_target_primers['id'])
    relative_target_pos = chosen_target_primers['TARGET_START'][0]
    chromosome = target.slice.chromosome

    ugs_left = create_ugs_fwd(chromosome,
                              target.slice.start_pos - (relative_target_pos - chosen_target_primers['LEFT_START'][0]),
                              target.slice.start_pos - (relative_target_pos - chosen_target_primers['LEFT_START'][0] -
                                                        chosen_target_primers['LEFT_START'][1])-1)

    ugs_right = create_ugs_rev(chromosome,
                               target.slice.start_pos + (chosen_target_primers['RIGHT_START'][0] - relative_target_pos -
                                                         chosen_target_primers['RIGHT_START'][1] + 1),
                               target.slice.start_pos + (chosen_target_primers['RIGHT_START'][0] -
                                                         relative_target_pos)
                               )

    new_slice = create_new_clice(chromosome,
                                 target.slice.start_pos -
                                 (relative_target_pos - chosen_target_primers['LEFT_START'][0]),
                                 target.slice.start_pos +
                                 (chosen_target_primers['RIGHT_START'][0] - relative_target_pos)
                                 )

    # checking to se that the UGS objects compatible with the Primer3 results
    try:
        if ugs_left.sequence.seq.decode('utf-8') != chosen_target_primers['LEFT'] or \
                        ugs_right.sequence.seq.decode('utf-8') != chosen_target_primers['RIGHT']:
            raise ValueError("The UGS seq does not equals primer3 results (target id: {})".format(target.name))

    except ValueError as ve:
        print(ve)

    te = TargetEnrichment.objects.create(
        chromosome=chromosome,
        left=ugs_left,
        right=ugs_right,
        planning_version=1,
    )
    ta = UMITargetedAmplicon.objects.create(
        left_ugs=ugs_left,
        right_ugs=ugs_right,
        umi_length=3,
        slice=new_slice,
    )

    return ta, te

