"""
This file demonstrates writing tests using the unittest module. These will pass
when you run "manage.py test".

Replace this with more appropriate tests for your application.
"""

from django.test import TestCase
from LinApp.models import Taxa, Assembly, TargetEnrichment, Target, Microsatellite

class queries_test(TestCase):
    def targets_enrichments_test(self):
        """
        Tests that 1 + 1 always equals 2.
        """
        taxa_object = Taxa.objects.get(taxonomy_id=10090)
        assembly = Assembly.objects.get(taxa=taxa_object, friendly_name='mm9')
        for enrichment in TargetEnrichment.objects.filter(left__assembly=assembly)\
                .prefetch_related('targets__microsatellite', 'physical_locations',
                                  'left__sequence', 'right__sequence'):
            for target in enrichment.targets.all():
                try:
                    mstarget = target.microsatellite
                    mstarget.repeat_type
                    mstarget.repeat_number
                    mstarget.repeat_unit
                except Microsatellite.DoesNotExist:
                    pass
                target.start_pos
                target.referencevalue.length
                enrichment.left.sequence.sequence
                enrichment.left.start_pos
                locations = enrichment.physical_locations.all()

def targets_tdv(request, taxa, assem):
    for enrichment in TargetEnrichment.objects.filter(left__assembly=assembly)\
            .prefetch_related('targets', 'targets__microsatellite', 'physical_locations', 'left', 'right', 'left__sequence', 'right__sequence'):
        for target in enrichment.targets.all():
            s += target.name + '\t'
            try:
                mstarget = target.microsatellite
                s += str(mstarget.repeat_type) + '\t'
                s += str(mstarget.repeat_number) + '\t'
                s += mstarget.repeat_unit + '\t'
            except Microsatellite.DoesNotExist:
                s += '\t\t\t'
            s += str(target.start_pos) + '\t'
            s += str(target.end_pos) + '\t'
            s += str(target.referencevalue.length) + '\t'
            s += enrichment.left.sequence.sequence + '\t'
            s += str(enrichment.left.start_pos) + '\t'
            s += str(enrichment.left.end_pos) + '\t'
            s += enrichment.right.sequence.sequence + '\t'
            s += str(enrichment.right.start_pos) + '\t'
            s += str(enrichment.right.end_pos) + '\t'
            locations = enrichment.physical_locations.all()
            if len(locations) == 1:
                s += locations[0].plate.name + '\t'
                s += locations[0].well + '\t'
            else:
                s += '\t\t'
            s += str(enrichment.passed_validation) + '\t'
            s += str(enrichment.validation_failure_id) + '\t'
            s += str(enrichment.validation_date) + '\t'
            s += '\r\n'
    return HttpResponse(s, content_type="text/plain")