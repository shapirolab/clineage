from misc.models import Taxa
from genomes.models import Assembly, Chromosome, DNASlice
from targeted_enrichment.planning.models import UGSPlus, UGSMinus, TargetEnrichment, Microsatellite
from targeted_enrichment.reagents.models import PCR1PrimerPairTER
from targeted_enrichment.amplicons.models import PCR1Amplicon
from primers.parts.models import DNABarcode1, DNABarcode2, \
    IlluminaReadingAdaptor1, IlluminaReadingAdaptor2, \
    IlluminaReadingAdaptor1Cuts, IlluminaReadingAdaptor2Cuts
from primers.synthesis.models import PCR1PlusPrimer, PCR1MinusPrimer

human_taxa, c = Taxa.objects.get_or_create(
    name="Homo sapiens",
    taxonomy_id=9606,
    rank='species',
    friendly_name="Human",
)

hg19_assembly, c = Assembly.objects.get_or_create(
    taxa=human_taxa,
    name='February 2009 Homo sapiens (GRCh37)',
    friendly_name='hg19',
)

hg19_chromosome, c = Chromosome.objects.get_or_create(
    assembly=hg19_assembly,
    name='X',
    sequence_length=155270560,
    cyclic=False,
)

slice_28727_left, c = DNASlice.objects.get_or_create(
    chromosome=hg19_chromosome,
    start_pos=81316094,
    end_pos=81316116,
)

slice_28727_target_a, c = DNASlice.objects.get_or_create(
    chromosome=hg19_chromosome,
    start_pos=81316131,
    end_pos=81316199,
)

slice_28727_target_b, c = DNASlice.objects.get_or_create(
    chromosome=hg19_chromosome,
    start_pos=81316201,
    end_pos=81316236,
)

slice_28727_right, c = DNASlice.objects.get_or_create(
    chromosome=hg19_chromosome,
    start_pos=81316243,
    end_pos=81316265,
)

slice_28734_left, c = DNASlice.objects.get_or_create(
    chromosome=hg19_chromosome,
    start_pos=54384674,
    end_pos=54384696,
)

slice_28734_target_a, c = DNASlice.objects.get_or_create(
    chromosome=hg19_chromosome,
    start_pos=54384788,
    end_pos=54384805,
)

slice_28734_right, c = DNASlice.objects.get_or_create(
    chromosome=hg19_chromosome,
    start_pos=54384807,
    end_pos=54384829,
)

ugs_28727_left, c = UGSPlus.objects.get_or_create(
    slice=slice_28727_left,
)

ugs_28727_right, c = UGSMinus.objects.get_or_create(
    slice=slice_28727_right,
)

ms_28727_a, c = Microsatellite.objects.get_or_create(
    id=1,
    name='X_81316131_81316199',
    slice=slice_28727_target_a,
    repeat_unit_len=3,
    repeat_unit_type='AAG',
    repeat_number=23,
    repeat_unit_ref_seq='TCT',
)

ms_28727_b, c = Microsatellite.objects.get_or_create(
    id=2,
    name='X_81316201_81316236',
    slice=slice_28727_target_b,
    repeat_unit_len=3,
    repeat_unit_type='AGC',
    repeat_number=12,
    repeat_unit_ref_seq='CTG',
)

te_28727, c = TargetEnrichment.objects.get_or_create(
    chromosome=hg19_chromosome,
    left=ugs_28727_left,
    right=ugs_28727_right,
    planning_version=1,
)
te_28727.targets = [ms_28727_a, ms_28727_b]
te_28727.save()

ugs_28734_left, c = UGSPlus.objects.get_or_create(
    slice=slice_28734_left,
)

ugs_28734_right, c = UGSMinus.objects.get_or_create(
    slice=slice_28734_right,
)

ms_28734_a, c = Microsatellite.objects.get_or_create(
    id=3,
    name='X_54384788_54384805',
    slice=slice_28734_target_a,
    repeat_unit_len=3,
    repeat_unit_type='AAG',
    repeat_number=6,
    repeat_unit_ref_seq='AGA',
)

dnabarcode1, c = DNABarcode1.objects.get_or_create(
    name='D710',
    _sequence='TCCGCGAA'
)

dnabarcode2, c = DNABarcode2.objects.get_or_create(
    name='D508',
    _sequence='GTACTGAC'
)

dnabarcode1_a, c = DNABarcode1.objects.get_or_create(
    name='D718',
    _sequence='TGGGAGCC'
)

dnabarcode2_a, c = DNABarcode2.objects.get_or_create(
    name='D502',
    _sequence='ATAGAGGC'
)

illuminareadingadaptor1, c = IlluminaReadingAdaptor1.objects.get_or_create(
    name='Illumina Standard Reading Adaptor1',
    _sequence='ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
)

illuminareadingadaptor2, c = IlluminaReadingAdaptor2.objects.get_or_create(
    name='Illumina Standard Reading Adaptor2',
    _sequence='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
)

illuminareadingadaptor1cuts, c = IlluminaReadingAdaptor1Cuts.objects.get_or_create(
    ira=illuminareadingadaptor1,
    overlap_start=11,
    overlap_end=33,
)

illuminareadingadaptor2cuts, c = IlluminaReadingAdaptor2Cuts.objects.get_or_create(
    ira=illuminareadingadaptor2,
    overlap_start=12,
    overlap_end=34,
)

primer_28727_left, c = PCR1PlusPrimer.objects.get_or_create(
    name='chrX_81316201_81316236_ACG_fwd',
    ugs=ugs_28727_left,
    irac=illuminareadingadaptor1cuts,
)

primer_28727_right, c = PCR1MinusPrimer.objects.get_or_create(
    name='Seq19272_rev',
    ugs=ugs_28727_right,
    irac=illuminareadingadaptor2cuts,
)

primer_28734_left, c = PCR1PlusPrimer.objects.get_or_create(
    name='Seq20808_fwd',
    ugs=ugs_28734_left,
    irac=illuminareadingadaptor1cuts,
)

primer_28734_right, c = PCR1MinusPrimer.objects.get_or_create(
    name='Seq20808_rev',
    ugs=ugs_28734_right,
    irac=illuminareadingadaptor2cuts,
)

te_28734, c = TargetEnrichment.objects.get_or_create(
    chromosome=hg19_chromosome,
    left=ugs_28734_left,
    right=ugs_28734_right,
    planning_version=1,
)
te_28734.targets = [ms_28734_a]
te_28734.save()

ter_28727, c = PCR1PrimerPairTER.objects.get_or_create(
    te=te_28727,
    passed_validation=False,
    left_primer=primer_28727_left,
    right_primer=primer_28727_right,
)

ter_28734, c = PCR1PrimerPairTER.objects.get_or_create(
    te=te_28734,
    passed_validation=True,
    left_primer=primer_28734_left,
    right_primer=primer_28734_right,
)

try:
    PCR1Amplicon.objects.get(ter=ter_28727)
except PCR1Amplicon.DoesNotExist:
    pu_28727 = PCR1Amplicon(
        id=1,
        ter=ter_28727,
    )
    pu_28727.infer_slice()
    pu_28727.save()

try:
    PCR1Amplicon.objects.get(ter=ter_28734)
except PCR1Amplicon.DoesNotExist:
    pu_28734 = PCR1Amplicon(
        id=2,
        ter=ter_28734,
    )
    pu_28734.infer_slice()
    pu_28734.save()