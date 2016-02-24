from plumbum import local

BOWTIE2BUILD_PATH = "/net/mraid11/export/data/dcsoft/home/Adam/Software/bowtie2-2.2.2/bowtie2-build"
bowtie2build = local[BOWTIE2BUILD_PATH]

# bowtie2-build \
#     ${fastq_folder}/Cluster/processing/${fastqName}/merged/${fastq_R1}_final_merged.fa \
#     ${fastq_folder}/Cluster/processing/${fastqName}/merged/${fastq_R1}_final_merged



BOWTIE2_PATH = "/net/mraid11/export/data/dcsoft/home/Adam/Software/bowtie2-2.2.2/bowtie2"
bowtie2 = local[BOWTIE2_PATH]
bowtie2_with_defaults = bowtie2["-p", "24",
                                "-a",
                                "--very-sensitive",]

# bowtie2 \
#     -p 24 \
#     -a \
#     --very-sensitive \
#     -x ${fastq_folder}/Cluster/processing/${fastqName}/merged/${fastq_R1}_final_merged \
#     -f ${TargetsFolder}/${TargetsName}_Primers.fa \
#         > ${fastq_folder}/Cluster/processing/${fastqName}/merged/${fastq_R1}_final_merged_Primers.sam

# /net/mraid11/export/data/dcsoft/home/Adam/Software/bowtie2-2.2.2/bowtie2 \
#     -a --very-sensitive \
#     -x index_1 \
#     -f /net/mraid11/export/dcstor/Ofir/ngs_fixtures/28727_and_28734_Primers.fa \
#     -S test2.sam