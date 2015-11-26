from linapp.models import *
from DBUtils import DBUtils
from django.db.models import Q
from django.contrib import admin
from mptt.admin import MPTTModelAdmin

class ExperimentAdmin(admin.ModelAdmin):

    def get_queryset(self, request):
        qs = super(ExperimentAdmin, self).get_queryset(request)
        if request.user.is_superuser:
            return qs
        else:
            return qs.filter(Q(is_public=True)|Q(users = request.user))

    def has_add_permission(self, request):
        return True #Everyone can add an experiment.

    def has_change_permission(self, request, obj=None):
        if not obj:
            return True # So they can see the change list page
        if request.user.is_superuser or DBUtils.getrole(request.user, obj).write:
            return True
        else:
            return False

    def has_delete_permission(self, request, obj=None):
        if request.user.is_superuser or DBUtils.getrole(request.user, obj).delete:
            return True
        else:
            return False

# #### Users/Roles Management
# admin.site.register(LineageRole)
#
# #### Experiment Management
# admin.site.register(FileContext)
#
# #### Types and descriptors
# admin.site.register(ProtocolType)
# admin.site.register(SampleComposition)#e.g. single cell or bulk
# admin.site.register(CellContentType)
# admin.site.register(TargetType)
# admin.site.register(TargetVariantType)
# admin.site.register(SampleStatus)
#
# #### Generic Biological Objects
# admin.site.register(Taxa)
# admin.site.register(Organ)
# admin.site.register(Tissue)
# admin.site.register(Assembly)
#
# #### Algorithms description
# admin.site.register(AlgorithmType)
#
# #### Sequencing
# admin.site.register(MachineType)
#
# #### Optional storage mapping
# admin.site.register(StorageType)
# admin.site.register(StorageBox)
# admin.site.register(PlateStorage)
# admin.site.register(PlateContext)
# admin.site.register(PlatePlastica)
# admin.site.register(PlateType)
#
# ####### Temporarily registered in admin interface for dev purposes ####### TODO: Remove those
# admin.site.register(Experiment, ExperimentAdmin)
# admin.site.register(UserProfile)
# admin.site.register(ExperimentUser)
# admin.site.register(Extraction)
# admin.site.register(SamplingEvent)
# admin.site.register(Coordinates)
# admin.site.register(CellContent)
# admin.site.register(Cell)
# admin.site.register(Sequencing)
# admin.site.register(Location)
# admin.site.register(Individual)
# admin.site.register(Panel)
# admin.site.register(GeneticBackground)
# admin.site.register(Protocol)
# admin.site.register(ExperimentFile)
# admin.site.register(Target)
# admin.site.register(Algorithm)
# admin.site.register(AlgorithmParameter)
# admin.site.register(AlgorithmRun)
# admin.site.register(AlgorithmRunParameters)
# admin.site.register(Machine)
# admin.site.register(CellTreeNode, MPTTModelAdmin)
# admin.site.register(RawData)
# admin.site.register(CorrectedRawData)


### Users/Roles Management
# FIXME admin.site.register(LineageRole)
# FIXME admin.site.register(UserProfile)

### Experiment Management
# FIXME admin.site.register(Experiment)
# FIXME admin.site.register(ExperimentLog)
# FIXME admin.site.register(FileContext)
# FIXME admin.site.register(ExperimentFile)
# FIXME admin.site.register(ExperimentUser)

### Types and descriptors
# FIXME admin.site.register(ProtocolType)
# FIXME admin.site.register(SampleComposition)#e.g. single cell or bulk
# FIXME admin.site.register(CellContentType)
# FIXME admin.site.register(TargetType)
# FIXME admin.site.register(TargetVariantType) #TODO: This might be useless. Decide.admin.site.register(SampleStatus)
# FIXME admin.site.register(TargetEnrichmentFailureType)

### Generic Biological Objects
# FIXME admin.site.register(Taxa)
# FIXME admin.site.register(GeneticBackground)
# FIXME admin.site.register(Organ)
# FIXME admin.site.register(Tissue)
# FIXME admin.site.register(Assembly)

#admin.site.register(Kit)
# FIXME admin.site.register(Protocol)
# FIXME admin.site.register(Sequence)
# FIXME admin.site.register(Target)#Target is a locus on a reference genome.
# FIXME admin.site.register(Primer)
# FIXME admin.site.register(Microsatellite)
# FIXME admin.site.register(TargetEnrichmentType)
# FIXME admin.site.register(TargetEnrichment)
# FIXME admin.site.register(Coordinates)

#XYZ coordinates of laser capture.
# FIXME admin.site.register(FACSMarker)
# FIXME admin.site.register(Panel)#collection of targets
# FIXME admin.site.register(Location)  # Freetext location

### Algorithms description
# FIXME admin.site.register(AlgorithmType)#e.g. raw data to consensus data
# FIXME admin.site.register(Algorithm)
# FIXME admin.site.register(AlgorithmParameter)
# FIXME admin.site.register(AlgorithmRun)
# FIXME admin.site.register(AlgorithmRunParameters)

### Sampling Hierarchy
# FIXME admin.site.register(Individual)
# FIXME admin.site.register(ExtractionEvent)
# FIXME admin.site.register(Extraction)
# FIXME admin.site.register(SamplingEvent)
# FIXME admin.site.register(FACS)
# FIXME admin.site.register(LaserCapture)
# FIXME admin.site.register(CellSelector)
# FIXME admin.site.register(Cell)
# FIXME admin.site.register(CellContent)  # aka DNA

### Sequencing
# FIXME admin.site.register(MachineType)
# FIXME admin.site.register(Machine)
# FIXME admin.site.register(Sequencing)

### Sequencing Data classes
# FIXME admin.site.register(RawData)
# FIXME admin.site.register(CorrectedRawData)

### Targets Hierarchy
# FIXME admin.site.register(FailedTargetValue)
# FIXME admin.site.register(SequenceDistribution)
# FIXME admin.site.register(TargetAnalysis)
# FIXME admin.site.register(TargetVariant)
# FIXME admin.site.register(GenSig)
# FIXME admin.site.register(DM)
# FIXME admin.site.register(CellTreeNode)

### Optional storage mapping
# FIXME admin.site.register(StorageType)
# FIXME admin.site.register(StorageBox)
# FIXME admin.site.register(PlateContext) #The plate's context in use. e.g. pcr
# FIXME admin.site.register(PlatePlastica) #The plate's physical form. e.g. deepwell square
# FIXME admin.site.register(PlateType)
# FIXME admin.site.register(Plate)
# FIXME admin.site.register(PlateStorage)
# FIXME admin.site.register(SampleLocation)
### Primers Multiplex
# FIXME admin.site.register(PrimersMultiplex)
