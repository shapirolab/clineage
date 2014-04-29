from LinApp.models import *
from DBUtils import DBUtils
from django.db.models import Q
from django.contrib import admin
from mptt.admin import MPTTModelAdmin

class ExperimentAdmin(admin.ModelAdmin):

    def queryset(self, request):
        qs = super(ExperimentAdmin, self).queryset(request)
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
admin.site.register(LineageRole)
admin.site.register(UserProfile)

### Experiment Management
admin.site.register(Experiment)
admin.site.register(ExperimentLog)
admin.site.register(FileContext)
admin.site.register(ExperimentFile)
admin.site.register(ExperimentUser)

### Types and descriptors
admin.site.register(ProtocolType)
admin.site.register(SampleComposition)#e.g. single cell or bulk
admin.site.register(CellContentType)
admin.site.register(TargetType)
admin.site.register(TargetVariantType) #TODO: This might be useless. Decide.admin.site.register(SampleStatus)
admin.site.register(TargetEnrichmentFailureType)

### Generic Biological Objects
admin.site.register(Taxa)
admin.site.register(GeneticBackground)
admin.site.register(Organ)
admin.site.register(Tissue)
admin.site.register(Assembly)

#admin.site.register(Kit)
admin.site.register(Protocol)
admin.site.register(Sequence)
admin.site.register(Target)#Target is a locus on a reference genome.
admin.site.register(Primer)
admin.site.register(Microsatellite)
admin.site.register(TargetEnrichmentType)
admin.site.register(TargetEnrichment)
admin.site.register(Coordinates)

#XYZ coordinates of laser capture.
admin.site.register(FACSMarker)
admin.site.register(Panel)#collection of targets
admin.site.register(Location)  # Freetext location

### Algorithms description
admin.site.register(AlgorithmType)#e.g. raw data to consensus data
admin.site.register(Algorithm)
admin.site.register(AlgorithmParameter)
admin.site.register(AlgorithmRun)
admin.site.register(AlgorithmRunParameters)

### Sampling Hierarchy
admin.site.register(Individual)
admin.site.register(ExtractionEvent)
admin.site.register(Extraction)
admin.site.register(SamplingEvent)
admin.site.register(FACS)
admin.site.register(LaserCapture)
admin.site.register(CellSelector)
admin.site.register(Cell)
admin.site.register(CellContent)  # aka DNA

### Sequencing
admin.site.register(MachineType)
admin.site.register(Machine)
admin.site.register(Sequencing)

### Sequencing Data classes
admin.site.register(RawData)
admin.site.register(CorrectedRawData)

### Targets Hierarchy
admin.site.register(FailedTargetValue)
admin.site.register(SequenceDistribution)
admin.site.register(TargetAnalysis)
admin.site.register(TargetVariant)
admin.site.register(GenSig)
admin.site.register(DM)
admin.site.register(CellTreeNode)

### Optional storage mapping
admin.site.register(StorageType)
admin.site.register(StorageBox)
admin.site.register(PlateContext) #The plate's context in use. e.g. pcr
admin.site.register(PlatePlastica) #The plate's physical form. e.g. deepwell square
admin.site.register(PlateType)
admin.site.register(Plate)
admin.site.register(PlateStorage)
admin.site.register(SampleLocation)
### Primers Multiplex
admin.site.register(PrimersMultiplex)