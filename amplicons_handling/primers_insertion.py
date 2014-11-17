__author__ = 'ofirr'

import re
import hashlib
import sys
from django.core.management import setup_environ
sys.path.append('/home/mint/clineage')
from clineage import settings
setup_environ(settings)
from linapp.models import Target, TargetEnrichment, Primer, TargetType, Sequence
from utils.SequenceManipulations import complement


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


class PrimerLocationError(Exception):
    def __init__(self, msg=None):
        if msg is None:
            msg = 'Cannot locate primer sequence within margins from target.'
        super(PrimerLocationError, self).__init__(msg)


class AmpliconCollisionError(Exception):
    def __init__(self, msg=None):
        if msg is None:
            msg = 'Exising overlapping amplicon.'
        super(AmpliconCollisionError, self).__init__(msg)


def check_primers(target, primer_left_sequence, primer_right_sequence, target_enrichment_type, margins=120):
    """
    Search given primer sequences surrounding a Target object in the reference genome
    of the target. Checks for existing overlapping amplicons.
    """
    try:
        pf_s, pf_e = target.chromosome.locate(target.start_pos,
                                              target.end_pos,
                                              primer_left_sequence,
                                              padding=margins)
        pr_s, pr_e = target.chromosome.locate(target.start_pos,
                                              target.end_pos,
                                              complement(primer_right_sequence)[::-1],
                                              padding=margins)
    except ValueError:
        raise PrimerLocationError
    colliding_te = TargetEnrichment.objects.filter(type=target_enrichment_type)\
                                      .filter(chromosome=target.chromosome)\
                                      .filter(left__start_pos__lte=pr_e)\
                                      .filter(right__end_pos__gte=pf_s)
    if colliding_te:
	print colliding_te
        raise AmpliconCollisionError
    return (pf_s, pf_e), (pr_s, pr_e)



def create_primers_in_db(chosen_target_primers, target_enrichment_type, primer_type=TargetType.objects.get(name='Flank'), pf_tail=None, pr_tail=None, margins=200):
    colliding_amplicons = []
    create_primer_pairs = []
    te_list = []
    for target_id in chosen_target_primers:
	print target_id
        target = Target.objects.get(pk=target_id)
        primer_left_sequence = chosen_target_primers[target_id]['LEFT']
        primer_right_sequence = chosen_target_primers[target_id]['RIGHT']
        try:
            primers_indexes_tuple = check_primers(target, primer_left_sequence, primer_right_sequence, target_enrichment_type, margins=margins)
        except PrimerLocationError:
            print 'Unresolved primers for target {}, pf:{}, pr:{}'.format(target.id, primer_left_sequence, primer_right_sequence)
            continue
        except AmpliconCollisionError:
	    print 'Colliding primers for target {}, pf:{}, pr:{}'.format(target.id, primer_left_sequence, primer_right_sequence)
            colliding_amplicons.append(target)
            continue
        left_primer_indexes, right_primer_indexes = primers_indexes_tuple
        pf_s, pf_e = left_primer_indexes
        pr_s, pr_e = right_primer_indexes
        pf_refseq = get_or_create_sequence(primer_left_sequence)
        pr_refseq = get_or_create_sequence(primer_right_sequence)
        if pf_tail and pr_tail:
            pf_seq = get_or_create_sequence(pf_tail.tail+primer_left_sequence)
            pr_seq = get_or_create_sequence(pr_tail.tail+complement(primer_right_sequence)[::-1])
        else:
            pf_seq = get_or_create_sequence(primer_left_sequence)
            pr_seq = get_or_create_sequence(complement(primer_right_sequence)[::-1])
        TargetType.objects.get(name='Flank')
        primer_fwd, created_fw = Primer.objects.get_or_create(start_pos=pf_s,
                                        end_pos=pf_e,
                                        defaults={'name': target.name + '_fwd',
                                                  'type': primer_type,
                                                  'chromosome': target.chromosome,
                                                  'referencevalue': pf_refseq,
                                                  'strand': '+',
                                                  'sequence': pf_seq,
                                                  'tail': pf_tail,
                                                  })
	print "Primer fw {} INFO: {}".format(primer_fwd, created_fw)
        primer_rev, created_rv = Primer.objects.get_or_create(start_pos=pr_s,
                                        end_pos=pr_e,
                                        defaults={'name': target.name + '_rev',
                                                  'type': primer_type,
                                                  'chromosome': target.chromosome,
                                                  'referencevalue': pr_refseq,
                                                  'strand': '-',
                                                  'sequence': pr_seq,
                                                  'tail': pr_tail,
                                                  })
	print "Primer rev {} INFO: {}".format(primer_rev, created_rv)
        te_made, created = TargetEnrichment.objects.get_or_create(
				    chromosome=target.chromosome,
				    left=primer_fwd,
				    right=primer_rev,
				    defaults={
				    'type': target_enrichment_type,
				    'amplicon': target.chromosome.getdna(primer_fwd.start_pos, primer_rev.end_pos),
				    'passed_validation': None,
				    'validation_failure': None,
				    'validation_date': None,
				    'comment': None,
				    })
	te_list.append(te_made)
    return te_list
