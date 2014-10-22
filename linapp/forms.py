from dojango.forms import *
#from django import forms
from linapp.models import *
from linapp.widgets import *
#### -------------------------------------------------------------------------------------
#### Users/Roles Management
#### -------------------------------------------------------------------------------------
class UserProfileForm(ModelForm):
    class Meta:
        model = UserProfile
#    user = models.OneToOneField(User)
#    institute = models.CharField(max_length=50)
#    comment = models.TextField()


#### -------------------------------------------------------------------------------------
#### Experiment Management
#### -------------------------------------------------------------------------------------
class ExperimentForm(ModelForm):
    class Meta:
        model = Experiment
#    users = models.ManyToManyField(User, related_name='experiments', through='ExperimentUser')
#    name = models.CharField(max_length=50)
#    created_date = models.DateField(auto_now_add=True) #Automatically set the field to now when the object is first created.
#    description = models.TextField()
#    is_public = models.BooleanField(default=False)

class ExperimentFileForm(ModelForm):
    class Meta:
        model = ExperimentFile
#    title = models.CharField(max_length=50)
#    experiment = models.ForeignKey(Experiment)
#    file_name = models.CharField(max_length=50)
#    file = models.FileField(upload_to=path)
#    context = models.ForeignKey(FileContext)
#    upload_date = models.DateField()
#    user = models.ForeignKey(User)
#    description = models.TextField()

class ExperimentLogForm(ModelForm):
    class Meta:
        model =ExperimentLog
#    user = models.ForeignKey(User)
#    experiment = models.ForeignKey(Experiment, related_name='comments')
#    date = models.DateTimeField(auto_now=True)
#    comment = models.TextField()

class ExperimentUserForm(ModelForm):
    class Meta:
        model =ExperimentUser
        fields = ('user', 'role')
#    user = models.ForeignKey(User)
#    experiment = models.ForeignKey(Experiment)
#    role = models.ForeignKey(LineageRole)


#### -------------------------------------------------------------------------------------
#### Generic Biological Objects
#### -------------------------------------------------------------------------------------
class GeneticBackgroundForm(ModelForm):
    class Meta:
        model =GeneticBackground
#    name = models.CharField(max_length=50)

class ProtocolForm(ModelForm):
    class Meta:
        model = Protocol
#    name = models.CharField(max_length=100)
#    abstract = models.TextField()
#    fulldescription = models.TextField()
#    kit = models.CharField(max_length=100)
#    type = models.ForeignKey(ProtocolType)

class TargetForm(ModelForm):
    class Meta:
        model = Target
#    name = models.CharField(max_length=50)
#    type = models.ForeignKey(TargetType) #Microsatellite / SNP / etc...
#    assembly = models.ForeignKey(Assembly)
#    chromosome = models.CharField(max_length=50)
#    start_pos = models.IntegerField()
#    end_pos = models.IntegerField()
#    referencevalue = models.TextField()

class CoordinatesForm(ModelForm):
    class Meta:
        model = Coordinates
#    x = models.DecimalField(max_digits=10, decimal_places=4)
#    y = models.DecimalField(max_digits=10, decimal_places=4)
#    z = models.DecimalField(max_digits=10, decimal_places=4)

class PanelForm(ModelForm):
    class Meta:
        model = Panel
#    name = models.CharField(max_length=50)
#    targets = models.ManyToManyField(Target, related_name='panels')

class LocationForm(ModelForm):
    class Meta:
        model = Location
#    name = models.CharField(max_length=50)


#### -------------------------------------------------------------------------------------
#### Algorithms description
#### -------------------------------------------------------------------------------------
class AlgorithmForm(ModelForm):
    class Meta:
        model = Algorithm
#    name = models.CharField(max_length=50)
#    type = models.ForeignKey(AlgorithmType)
#    version = models.CharField(max_length=50)
#    developers = models.ManyToManyField(User)

class AlgorithmParameterForm(ModelForm):
    class Meta:
        model = AlgorithmParameter
#    algorithm = models.ManyToManyField(Algorithm, related_name='parameters')
#    name = models.CharField(max_length=50)

class AlgorithmRunForm(ModelForm):
    class Meta:
        model = AlgorithmRun
#    algorithm = models.ForeignKey(Algorithm, related_name='runs')
#    runname = models.CharField(max_length=50, null=True)
#    parameters = models.ManyToManyField(AlgorithmParameter, through='AlgorithmRunParameters')
#    user = models.ForeignKey(User)
#    timestamp = models.DateTimeField()
#    ExtraFiles =models.ManyToManyField(ExperimentFile)

class AlgorithmRunParametersForm(ModelForm):
    class Meta:
        model = AlgorithmRunParameters
#    run = models.ForeignKey(AlgorithmRun)
#    parameter = models.ForeignKey(AlgorithmParameter)
#    value = models.CharField(max_length=50)


#### -------------------------------------------------------------------------------------
#### Sampling Hierarchy
#### -------------------------------------------------------------------------------------
class IndividualForm(ModelForm):
    class Meta:
        model = Individual
#    GENDER = (('M', 'Male'),('F', 'Female'),)
#    taxa = models.ForeignKey(Taxa)
#    name = models.CharField(max_length=50)
#    sex = models.CharField(max_length=1, choices=GENDER)
#    born = models.DateTimeField()
#    comment = models.TextField()
#    background = models.ForeignKey(GeneticBackground)
#    location = models.ForeignKey(Location,null=True)

class ExtractionForm(ModelForm):
    class Meta:
        model = Extraction
#    individual = models.ForeignKey(Individual)
#    name = models.CharField(max_length=50)
#    user = models.ForeignKey(User)
#    date = models.DateTimeField()
#    location = models.ForeignKey(Location,null=True)
#    comment = models.TextField()

class SamplingEventForm(ModelForm):
    class Meta:
        model = SamplingEvent
#    name = models.CharField(max_length=50)
#    organ = models.ForeignKey(Organ)
#    extraction = models.ForeignKey(Extraction)
#    coordinates = models.ForeignKey(Coordinates)
#    tissue = models.ForeignKey(Tissue)
#    user = models.ForeignKey(User)
#    comment = models.TextField()

class CellForm(ModelForm):
    class Meta:
        model = Cell
#    sampling = models.ForeignKey(SamplingEvent)
#    name = models.CharField(max_length=50)
#    experiment = models.ManyToManyField(Experiment, related_name='cells')
#    composition = models.ForeignKey(SampleComposition)#single cell or bulk
#    status = models.ForeignKey(SampleStatus)
#    comment = models.TextField()

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
    # def __init__(self, *args, **kwargs):
    #     # call constructor to set up the fields. If you don't do this
    #     # first you can't modify fields.
    #     super(PlateInputForm, self).__init__(*args, **kwargs)
    #     try:
    #         revko = StorageType.objects.get(name='REVCO')
    #         freezer = StorageType.objects.get(name='Freezer')
    #         if self.instance.storage_box.storage_type == revko:
    #             inner_location = RevcoWidget()
    #         elif self.instance.storage_box.storage_type == freezer:
    #             inner_location = FreezerWidget()
    #         else:
    #             inner_location = HiddenInput()
    #     except AttributeError:
    #         # unbound form
    #         inner_location = HiddenInput()

    # cell input fields
    individual = ModelChoiceField(queryset=Individual.objects.all())
    sampling = ModelChoiceField(queryset=SamplingEvent.objects.all(), required=False)
    cells_name_prefix = CharField(max_length=50, required=False)
    comment = CharField(widget=Textarea, required=False)

    # plate input fields
    existing_plate = ModelChoiceField(queryset=Plate.objects.all(), required=False)
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


class CellContentForm(ModelForm):
    class Meta:
        model = CellContent
#    parent = models.ForeignKey('CellContent', null=True, blank=True)
#    cell = models.ForeignKey(Cell)
#    panel = models.ForeignKey(Panel, null=True, blank=True)
#    type = models.ForeignKey(CellContentType)
#    name = models.CharField(max_length=50)
#    protocol = models.ForeignKey(Protocol)
#    user = models.ForeignKey(User)
#    comment = models.TextField()


#### -------------------------------------------------------------------------------------
#### Sequencing
#### -------------------------------------------------------------------------------------
class MachineForm(ModelForm):
    class Meta:
        model = Machine
#    machineid = models.CharField(max_length=50)
#    type = models.ForeignKey(MachineType)

class SequencingForm(ModelForm):
    class Meta:
        model = Sequencing
#    sample = models.ForeignKey(CellContent)
#    data = models.ForeignKey('RawData', related_name='sequencing_event') #we add '' brackets in order to resolve precedence issue
#    machine = models.ForeignKey(Machine)
#    protocol = models.ForeignKey(Protocol)
#    user = models.ForeignKey(User)
#    date = models.DateField()


#### -------------------------------------------------------------------------------------
#### Sequencing Data classes
#### -------------------------------------------------------------------------------------
class RawDataForm(ModelForm):
    class Meta:
        model = RawData
#    sequencing = models.ForeignKey(Sequencing)
#    file = models.ForeignKey(ExperimentFile)

class CorrectedRawDataForm(ModelForm):
    class Meta:
        model = CorrectedRawData
#    rawdata = models.ForeignKey(RawData)
#    file = models.ForeignKey(ExperimentFile)
#    creation = models.ForeignKey(AlgorithmRun, related_name = 'crd')


#### -------------------------------------------------------------------------------------
#### Targets Hierarchy
#### -------------------------------------------------------------------------------------
class FailedTargetValueForm(ModelForm):
    class Meta:
        model = FailedTargetValue
#    comment = models.TextField()

class SequenceDistributionForm(ModelForm):
    class Meta:
        model = SequenceDistribution
#    target = models.ManyToManyField(Target, through = 'TargetAnalysis')
#    value = models.TextField()
#    count = models.IntegerField()
#    failed = models.ForeignKey(FailedTargetValue, null=True)

class TargetAnalysisForm(ModelForm):
    class Meta:
        model = TargetAnalysis
#    target = models.ForeignKey(Target)
#    distribution = models.ForeignKey(SequenceDistribution)
#    creation = models.ForeignKey(AlgorithmRun, related_name='targetsdistributions')
#    sequencingdata = models.ForeignKey(CorrectedRawData)

class TargetVariantForm(ModelForm):
    class Meta:
        model = TargetVariant
#    ts = models.ForeignKey(TrueSequence)
#    creation = models.ForeignKey(AlgorithmRun,related_name='targetsvariants')
#    type = models.ForeignKey(TargetVariantType)
#    value = models.DecimalField(max_digits=10, decimal_places=4)

class GenSigForm(ModelForm):
    class Meta:
        model = GenSig
#    variants = models.ManyToManyField(TargetVariant)
#    creation = models.ForeignKey(AlgorithmRun,related_name='gensigs')
#    value = models.TextField() #vector

class DMForm(ModelForm):
    class Meta:
        model = DM
#    cell1 = models.ForeignKey(GenSig, related_name='+')
#    cell2 = models.ForeignKey(GenSig, related_name='+')
#    distance = models.DecimalField(max_digits=10, decimal_places=4)
#    creation = models.ForeignKey(AlgorithmRun, related_name='dms')

class CellTreeNodeForm(ModelForm): #TODO: This might be MPTTForm
    class Meta:
        model = CellTreeNode
#    cell = models.ForeignKey(Cell)
#    parent = TreeForeignKey('self', null=True, blank = True, related_name='children')
#    distance = models.DecimalField(max_digits=10, decimal_places=4)
#    def __unicode__(self):
#        return self.cell.name
#    #    class MPTTMeta:
##        order_insertion_by=['name']


#### -------------------------------------------------------------------------------------
#### Optional storage mapping
#### -------------------------------------------------------------------------------------
class PlateForm(ModelForm):
    class Meta:
        model = Plate
#    code = models.AutoField(primary_key=True)
#    type = models.ForeignKey(PlateType)
#    name = models.CharField(max_length=200, blank=True)
#    barcode = models.CharField(max_length=20, blank=True)
#    timestamp = models.DateField(null=True, blank=True)
#    state = models.CharField(max_length=20, blank=True)
#    lastusedwell = models.CharField(max_length=4)

class PlateStorageForm(ModelForm):
    class Meta:
        model = PlateStorage
#    storageBox = models.ForeignKey(StorageBox)
#    plate = models.ForeignKey(Plate)
#    notes = models.CharField(max_length=250, blank=True)

class SampleLocationForm(ModelForm):
    class Meta:
        model = SampleLocation
#    plate = models.ForeignKey(Plate)
#    well = models.CharField(max_length=3, blank=True)
#    sample = models.ForeignKey(CellContent)
#    volume = models.DecimalField(null=True, max_digits=10, decimal_places=3, blank=True)
#    concentration = models.DecimalField(null=True, max_digits=10, decimal_places=5, blank=True)
#    timestamp = models.DateTimeField(auto_now=True)

# class PlateCellsForm(Form):