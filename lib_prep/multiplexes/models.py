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

    def __str__(self):
        return self.name


class PCR1Multiplex(TERMultiplex):
    ters = models.ManyToManyField(PCR1PrimerPairTERBase)
    #physical_locations = fields.GenericRelation(SampleLocation,
                               #content_type_field='content_type',
                               #object_id_field='object_id')
    #TODO: physical_locations(MPXPlate)


class PCR1Panel(models.Model):  # PCR1MultiplexCollection
    name = models.CharField(max_length=20)
    mpxs = models.ManyToManyField(PCR1Multiplex)
