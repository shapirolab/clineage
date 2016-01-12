
from django.db import models


### -------------------------------------------------------------------------------------
### Generic Biological Objects
### -------------------------------------------------------------------------------------


class Taxa(models.Model):
    name = models.CharField(max_length=50)
    taxonomy_id = models.IntegerField()
    rank = models.CharField(max_length=50)
    parent = models.IntegerField(null=True, blank=True)
    friendly_name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.name

### -------------------------------------------------------------------------------------

