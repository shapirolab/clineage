__author__ = 'ofirr'

from linapp.models import TargetEnrichment

def other_target_enrichment_withing_margins(target_enrichment, existing_target_enrichments, margins=10000):
    """Check if the new amplicon collides with one of the existing amplicons within given margins"""
    existing_target_enrichment_pks = [te.pk for te in existing_target_enrichments]
    existing_target_enrichments_query = TargetEnrichment.objects.filter(pk__in=existing_target_enrichment_pks)
    return existing_target_enrichments_query.filter(left__start_pos__lte=target_enrichment.right.end_pos+margins).\
        filter(right__end_pos__gte=target_enrichment.left.start_pos-margins)


def populate_conditioned_mpx(conditioned_mpx, target_enrichment, margins=10000):
    """Check if the new amplicon satisfies all mpx conditions and add it"""
    if conditioned_mpx['size'] <= len(conditioned_mpx['target_enrichments']):
        return False, conditioned_mpx
    if other_target_enrichment_withing_margins(target_enrichment, conditioned_mpx['target_enrichments'], margins=margins):
        return False, conditioned_mpx
    if target_enrichment in conditioned_mpx['allowed_target_enrichments']:
        conditioned_mpx['target_enrichments'].append(target_enrichment)
        return True, conditioned_mpx
    return False, conditioned_mpx


def populate_conditioned_mpxs_map(target_enrichment_queryset, conditioned_mpxs, margins=10000):
    """Iterate over available target enrichments and add them to the first multiplex group who's conditions are all satisfied"""
    misfits = []
    for target_enrichment in target_enrichment_queryset:
        added = False
        for conditioned_mpx in conditioned_mpxs.keys():
            added, conditioned_mpxs[conditioned_mpx] = populate_conditioned_mpx(conditioned_mpxs[conditioned_mpx], target_enrichment, margins=margins)
            if added:
                break
        if not added:
            misfits.append(target_enrichment)
        temp_sum = reduce(lambda x,y:x+y, [len(conditioned_mpxs[group]['target_enrichments']) for group in conditioned_mpxs.keys()])
        expected_sum = reduce(lambda x,y:x+y, [conditioned_mpxs[group]['size'] for group in conditioned_mpxs.keys()])
        if temp_sum == expected_sum:
            return misfits, conditioned_mpxs
        if temp_sum > expected_sum:
            return misfits, conditioned_mpxs  # This should never happen
    return misfits, conditioned_mpxs