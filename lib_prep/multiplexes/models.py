from django.contrib.contenttypes import fields
from django.db import models
from targeted_enrichment.planning.models import TargetEnrichment
from targeted_enrichment.reagents.models import PCR1PrimerPairTERBase
from wet_storage.models import SampleLocation

### -------------------------------------------------------------------------------------
### Primers Multiplex
### -------------------------------------------------------------------------------------


#TODO : notail?


class TERMultiplex(models.Model): # TODO: move to primers, m2m to TER.
    name = models.CharField(max_length=20)
    # ters = models.ManyToManyField(TargetedEnrichmentReagent)
    physical_locations = fields.GenericRelation(SampleLocation,
                               content_type_field='content_type',
                               object_id_field='object_id')

    class Meta:
        abstract = True

    def __unicode__(self):
        return self.name


class PCR1Multiplex(TERMultiplex):
    ters = models.ManyToManyField(PCR1PrimerPairTERBase)
    #physical_locations = fields.GenericRelation(SampleLocation,
                               #content_type_field='content_type',
                               #object_id_field='object_id')
    #TODO: physical_locations(MPXPlate)


class Panel(models.Model):#collection of targets
                                                # TODO: m2m pri_mux, well on the m2m table.
    name = models.CharField(max_length=50)
    targets = models.ManyToManyField(TargetEnrichment, related_name='panels')
    def __unicode__(self):
        return self.name
