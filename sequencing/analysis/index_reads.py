from plumbum import local

BOWTIE2BUILD_PATH = "/net/mraid11/export/data/dcsoft/home/Adam/Software/bowtie2-2.2.2/bowtie2-build"
bowtie2build = local[BOWTIE2BUILD_PATH]
