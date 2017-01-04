
import pytest

from tests.sequencing.analysis.adamiya.conftest import *
from tests.sequencing.analysis.adamiya.conftest import _chain_histogram_entry_reads, _chain_amplicon_reads


@pytest.mark.skipif(pytest.config.getoption("nomigrations"), reason="No migrations, no view.")
@pytest.mark.django_db
def test_trigger_histogram_entry_reads_trg(adam_histogram_entry_reads_files_d, _chain_histogram_entry_reads, requires_microsatellites, requires_none_genotypes):
    second_num_read = 7
    none_snp_genotype = SNPHistogramGenotype.objects.get(snp=None)
    snp_histogram_genotypes, c = SNPHistogramGenotypeSet.objects.get_or_create(
        **{fn: none_snp_genotype for fn in SNPHistogramGenotypeSet.genotype_field_names()})
    for (l_id, bc_id, inc, amp), f_d_d in adam_histogram_entry_reads_files_d.items():
        h = _chain_histogram_entry_reads[l_id, bc_id, inc, amp]
        hid = h.id
        total_num_reads = 0
        for msgs, f_d in f_d_d.items():
            msg_objs = set()
            for (microsatellite_id, repeat_number) in msgs:
                msg, c = MicrosatelliteHistogramGenotype.objects.get_or_create(
                    microsatellite_id=microsatellite_id,
                    repeat_number=repeat_number,
                )
                msg_objs.add(msg)
            ms_genotypes = MicrosatelliteHistogramGenotypeSet.get_for_msgs(msg_objs)
            f_d2 = {}
            for r in [R1, R2, RM]:
                f_d2[r] = get_unique_path("fastq")
                os.symlink(f_d[r], f_d2[r])
            her, c = HistogramEntryReads.objects.get_or_create(
                histogram=h,
                microsatellite_genotypes=ms_genotypes,
                snp_genotypes=snp_histogram_genotypes,
                defaults=dict(num_reads=f_d[NUM_READS],
                fastqm=f_d2[RM],
                fastq1=f_d2[R1],
                fastq2=f_d2[R2]
            ))
            assert c
            total_num_reads += f_d[NUM_READS]
            assert total_num_reads == Histogram.objects.get(id=hid).num_reads

            # remove objects
            os.unlink(her.fastqm)
            os.unlink(her.fastq1)
            os.unlink(her.fastq2)


            msg_objs = set()
            for (microsatellite_id, repeat_number) in msgs:
                msg, c = MicrosatelliteHistogramGenotype.objects.get_or_create(
                    microsatellite_id=microsatellite_id,
                    repeat_number=repeat_number+3, #+3 to make sure a new ms is created
                )
                msg_objs.add(msg)
            ms_genotypes = MicrosatelliteHistogramGenotypeSet.get_for_msgs(msg_objs)
            f_d2 = {}
            for r in [R1, R2, RM]:
                f_d2[r] = get_unique_path("fastq")
                os.symlink(f_d[r], f_d2[r])
            her, c = HistogramEntryReads.objects.get_or_create(
                histogram=h,
                microsatellite_genotypes=ms_genotypes,
                snp_genotypes=snp_histogram_genotypes,
                defaults=dict(num_reads=second_num_read,
                fastqm=f_d2[RM],
                fastq1=f_d2[R1],
                fastq2=f_d2[R2]
            ))
            assert c
            total_num_reads += second_num_read
            assert total_num_reads == Histogram.objects.get(id=hid).num_reads
            
            #remove objects
            os.unlink(her.fastqm)
            os.unlink(her.fastq1)
            os.unlink(her.fastq2)





