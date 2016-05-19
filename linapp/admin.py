
from django.contrib import admin

from wet_storage.models import StorageType, StorageBox, PlateContext, \
    PlatePlastica, PlateType, Plate, PlateStorage, SampleLocation
from linapp.models import LineageRole, UserProfile, ProtocolType, Protocol, \
    UserReport
from lib_prep.workflows.models import BarcodedContent, CellContentProtocol, \
    AmplifiedContent
from lib_prep.multiplexes.models import PCR1Multiplex, PCR1Panel
from misc.models import Taxa
from genomes.models import Assembly, Chromosome, DNASlice
from targeted_enrichment.reagents.models import TargetEnrichmentFailureType, \
    PCR1PrimerPairTER, PCR1WithCompanyTagPrimerPairTER, \
    PCR1PrimerPairTERDeprecated, TargetedNoTailPrimerPairTER
from targeted_enrichment.planning.models import UGSPlus, UGSMinus, Target, \
    TargetEnrichment, RestrictionEnzyme, RestrictionSite, Microsatellite, SNP
from sampling.models import GeneticBackground, Organ, Tissue, \
    SampleComposition, SampleStatus, Coordinates, FACSMarker, Location, \
    Individual, ExtractionEvent, Extraction, SamplingEvent, FACS, \
    LaserCapture, CellSelector, Cell
from primers.parts.models import IlluminaReadingAdaptor1, \
    IlluminaReadingAdaptor2, IlluminaReadingAdaptor1ForTail, \
    IlluminaReadingAdaptor1ForHead, IlluminaReadingAdaptor2ForTail, \
    IlluminaReadingAdaptor2ForHead, IlluminaFlowCellAdaptor1, \
    IlluminaFlowCellAdaptor2, DNABarcode1, DNABarcode2
from primers.synthesis.models import TargetedNoTailPlusPrimer, \
    TargetedNoTailMinusPrimer, PCR1PlusPrimer, PCR1MinusPrimer, \
    PCR1WithCompanyTagPlusPrimer, PCR1WithCompanyTagMinusPrimer, \
    PCR2PlusPrimer, PCR2MinusPrimer
from sequencing.runs.models import MachineType, Machine, NGSRun

admin.site.register(StorageType)
admin.site.register(StorageBox)
admin.site.register(PlateContext)
admin.site.register(PlatePlastica)
admin.site.register(PlateType)
admin.site.register(Plate)
admin.site.register(PlateStorage)
admin.site.register(SampleLocation)

admin.site.register(LineageRole)
admin.site.register(UserProfile)
admin.site.register(ProtocolType)
admin.site.register(Protocol)
admin.site.register(UserReport)

admin.site.register(BarcodedContent)
admin.site.register(CellContentProtocol)
admin.site.register(AmplifiedContent)

admin.site.register(PCR1Multiplex)
admin.site.register(PCR1Panel)

admin.site.register(Taxa)

admin.site.register(Assembly)
admin.site.register(Chromosome)
admin.site.register(DNASlice)

admin.site.register(TargetEnrichmentFailureType)
admin.site.register(PCR1PrimerPairTER)
admin.site.register(PCR1WithCompanyTagPrimerPairTER)
admin.site.register(PCR1PrimerPairTERDeprecated)
admin.site.register(TargetedNoTailPrimerPairTER)

admin.site.register(UGSPlus)
admin.site.register(UGSMinus)
admin.site.register(Target)
admin.site.register(TargetEnrichment)
admin.site.register(RestrictionEnzyme)
admin.site.register(RestrictionSite)
admin.site.register(Microsatellite)
admin.site.register(SNP)

admin.site.register(GeneticBackground)
admin.site.register(Organ)
admin.site.register(Tissue)
admin.site.register(SampleComposition)
admin.site.register(SampleStatus)
admin.site.register(Coordinates)
admin.site.register(FACSMarker)
admin.site.register(Location)
admin.site.register(Individual)
admin.site.register(ExtractionEvent)
admin.site.register(Extraction)
admin.site.register(SamplingEvent)
admin.site.register(FACS)
admin.site.register(LaserCapture)
admin.site.register(CellSelector)
admin.site.register(Cell)

admin.site.register(IlluminaReadingAdaptor1)
admin.site.register(IlluminaReadingAdaptor2)
admin.site.register(IlluminaReadingAdaptor1ForTail)
admin.site.register(IlluminaReadingAdaptor1ForHead)
admin.site.register(IlluminaReadingAdaptor2ForTail)
admin.site.register(IlluminaReadingAdaptor2ForHead)
admin.site.register(IlluminaFlowCellAdaptor1)
admin.site.register(IlluminaFlowCellAdaptor2)
admin.site.register(DNABarcode1)
admin.site.register(DNABarcode2)

admin.site.register(TargetedNoTailPlusPrimer)
admin.site.register(TargetedNoTailMinusPrimer)
admin.site.register(PCR1PlusPrimer)
admin.site.register(PCR1MinusPrimer)
admin.site.register(PCR1WithCompanyTagPlusPrimer)
admin.site.register(PCR1WithCompanyTagMinusPrimer)
admin.site.register(PCR2PlusPrimer)
admin.site.register(PCR2MinusPrimer)

admin.site.register(MachineType)
admin.site.register(Machine)
admin.site.register(NGSRun)
