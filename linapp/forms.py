from dojango.forms import *
#from django import forms
from django.contrib.auth.models import User
from sampling.models import Individual, SamplingEvent, SampleComposition, \
    SampleStatus
from targeted_enrichment.planning.models import Target
from lib_prep.multiplexes.models import Panel
from wet_storage.models import Plate, PlateType, StorageBox
from linapp.models import Protocol 
from linapp.widgets import *


class MultipleCellForm(forms.Form):
    copies = IntegerField(min_value=1)
    individual = ModelChoiceField(queryset=Individual.objects.all())
    sampling = ModelChoiceField(queryset=SamplingEvent.objects.all(), required=False)
    cells_name_prefix = CharField(max_length=50, required=False)
    # experiment = ModelMultipleChoiceField(queryset=Experiment.objects.all())
    composition = ModelChoiceField(queryset=SampleComposition.objects.all())  # single cell or bulk
    status = ModelChoiceField(queryset=SampleStatus.objects.all(), required=False)
    comment = CharField(widget=Textarea, required=False)


class PlateInputForm(forms.Form):

    # cell input fields
    individual = ModelChoiceField(queryset=Individual.objects.all())
    sampling = ModelChoiceField(queryset=SamplingEvent.objects.all(), required=False)
    cells_name_prefix = CharField(max_length=50, required=False)
    comment = CharField(widget=Textarea, required=False)

    # plate input fields
    existing_plate = ModelChoiceField(
        queryset=Plate.objects.exclude(name__startswith='MM_v2_P')\
                              .exclude(name__startswith='amplicon_MM_notails_plt')\
                              .exclude(name__startswith='Primers_')\
                              .exclude(name__startswith='United_')\
                              .exclude(name__startswith='MM_MPX_w/o_LB_P')\
                              .exclude(name__startswith='AAR')\
                              .exclude(name__startswith='hg19_Tails_plt'),
        required=False)
    plate_type = ModelChoiceField(queryset=PlateType.objects.all(), required=False)
    plate_name = CharField(max_length=200, required=False)
    samples_in_wells = CharField(widget=Textarea)

    #plate storage fields
    # storage_type = ModelChoiceField(queryset=StorageType.objects.all())
    storage_box = ModelChoiceField(queryset=StorageBox.objects.select_related('storage_type'))
    inner_location = CharField(max_length=100, required=False)
    notes = CharField(max_length=250, required=False)

    # logical fields
    user = ModelChoiceField(User.objects.all())
    timestamp = DateField()
    protocol = ModelChoiceField(Protocol.objects.all(), required=False)
