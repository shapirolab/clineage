import os
import csv
import seaborn as sns  # development
from collections import defaultdict

from django.utils.encoding import smart_text
from django.contrib.auth.models import User

from targeted_enrichment.planning.models import Microsatellite, SNP
from lib_prep.multiplexes.models import OM6Panel, PCR1Panel
from targeted_enrichment.reagents.models import PCR1PrimerPairTER, OM6PadlockTER

def query_panel(panel_name='all'):
    panel_list = list()
    if panel_name=='all':
        for pcr1panel in PCR1Panel.objects.all():
            panel_list.append(pcr1panel)
        for om6panel in OM6Panel.objects.all():
            panel_list.append(om6panel)
    else:
        for pcr1panel in PCR1Panel.objects.filter(name__contains=panel_name):
            panel_list.append(pcr1panel)
        for om6panel in OM6Panel.objects.filter(name__contains=panel_name):
            panel_list.append(om6panel)
        return panel_list
    return panel_list


def panel_ms_cell_table_values_db(panel_name='all'):
    panels=query_panel(panel_name)
    for panel in panels:
        if type(panel)== OM6Panel:
            mxs= panel.mixs.all()
        if type(panel)== PCR1Panel:
            mxs= panel.mpxs.all()
        for mx in mxs:
            for ter in mx.ters.select_subclasses():
                amplicon = ter.amplicon
                for slice_contains in amplicon.slice.contains.all():
                    for target in slice_contains.target_set.all():
                        snp = None
                        ms = None
                        try:
                            ms=target.microsatellite
                        except Microsatellite.DoesNotExist:
                            try:
                                snp=target.snp
                            except SNP.DoesNotExist:
                                continue
                        yield{
                            'MS ID': smart_text(ms.id if ms is not None else ''),
                            'SNP ID': smart_text(snp.id if snp is not None else ''),
                            'Target ID': smart_text(target.id if target is not None else ''),
                            'Old Adam TEpk': smart_text(ter.old_adam_te_pk),
                            'TER ID':smart_text(ter.id),
                            'Amplicon ID': smart_text(amplicon.id),
                            'Basic Unit size': smart_text(ms.repeat_unit_len if ms is not None else ''),
                            'Expected Number of repeats': smart_text(ms.repeat_number if ms is not None else ''),
                            'Basic Unit Type': smart_text(ms.repeat_unit_type if ms is not None else ''),
                            'SNP Sequence':smart_text(snp.slice.sequence.seq.decode('utf-8') if snp is not None else ''),
                            'SNP modified':smart_text(snp.modified if snp is not None else ''),
                            'Chromosome': smart_text(ms.slice.chromosome if ms is not None else snp.slice.chromosome),
                            'Length MS': smart_text(abs(ms.slice.end_pos - ms.slice.start_pos) if ms is not None else ''),
                            'Target sequence': smart_text(target.slice.sequence.seq.decode('utf-8')),
                            'Primer sequence - Fw': smart_text(ter.left_primer.sequence.seq.decode('utf-8') \
                                       if type(panel) == PCR1Panel else ter.padlock.left_ugs.sequence.seq.decode('utf-8')),
                            'Primer sequence -  Rev': smart_text(ter.right_primer.sequence.seq.decode('utf-8')\
                                       if type(panel) == PCR1Panel else ter.padlock.right_ugs.sequence.seq.decode('utf-8')),
                            'Target location on Chromosome - start': smart_text(ter.te.left.slice.start_pos),
                            'Target location on Chromosome - end': smart_text(ter.te.right.slice.end_pos),
                            'Fw primer location on Chromosome - start': smart_text(ter.left_primer.ugs.slice.start_pos \
                                        if type(panel)== PCR1Panel else ter.padlock.left_ugs.slice.start_pos),
                            'Fw primer location on Chromosome - end': smart_text(ter.left_primer.ugs.slice.end_pos \
                                        if type(panel) == PCR1Panel else ter.padlock.left_ugs.slice.end_pos),
                            'Rev primer location on Chromosome - start': smart_text(ter.right_primer.ugs.slice.start_pos \
                                        if type(panel) == PCR1Panel else ter.padlock.right_ugs.slice.start_pos),
                            'Rev primer location on Chromosome - end': smart_text(ter.right_primer.ugs.slice.end_pos \
                                        if type(panel) == PCR1Panel else ter.padlock.right_ugs.slice.end_pos),
                            'Amplicon location on Chromosome - start': smart_text(amplicon.slice.start_pos,),
                            'Amplicon location on Chromosome - end': smart_text(amplicon.slice.end_pos),
                            'Amplicon Sequence': smart_text(amplicon.slice.sequence.seq.decode('utf-8')),
                            'MS location on Chromosome - start': smart_text(ms.slice.start_pos if ms is not None else ''),
                            'MS location on Chromosome - end': smart_text(ms.slice.end_pos if ms is not None else ''),
                            'SNP location on the Chromosome':smart_text(snp.slice.start_pos if snp is not None else ''),
                            'Mpx groups names': smart_text(mx.name if type(panel)== PCR1Panel else ''),
                            'Mpx group size': smart_text(len(mx.ters.select_subclasses()) if type(panel)== PCR1Panel else ''),
                            'Mixs groups names':smart_text(mx.name if type(panel)== OM6Panel else ''),
                            'Mixs group size': smart_text(len(mx.ters.select_subclasses()) if type(panel)== OM6Panel else '')

                        }