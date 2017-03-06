# __author__ = 'veronika'

import subprocess
import os
import plumbum

from misc.utils import get_unique_path, unlink, unique_file_cm
from Bio.Sequencing.Applications import SamtoolsMpileupCommandline

samtools = plumbum.local["samtools"]
samtools_mpileup = samtools["mpileup"]
samtools_mpileup_with_defaults = samtools_mpileup["-uvA",
                                                  "-d", "10000000",
                                                  "-t", "DP"]


def bcf_tools(bcf_mpileup_file):

    bashCommand = 'bcftools view -vcgI {}'.format(bcf_mpileup_file)
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

    vcf_file = get_unique_path(ext='vcf')
    str_data = process.communicate()[0].decode("utf-8")
    with open(vcf_file, 'w') as f:
        f.write(str_data)

    os.unlink(bcf_mpileup_file)
    return vcf_file


def mpileup(sort_bam_file, ref_faidx, d=10000000, A=True, u=True, D=True):

    # tmp_bcf_mpileup_data = SamtoolsMpileupCommandline(input_file=[sort_bam_file], f=ref_faidx, d=d, A=A, u=u, D=D)

    # bcf_mpileup_file = unique_file_cm(ext='bcf')
    # tmp_bcf_file_stderr = unique_file_cm(ext='txt')
    # tmp_bcf_mpileup_data(stdout=bcf_mpileup_file, stderr=tmp_bcf_file_stderr)

    # bashCommand = 'samtools mpileup -uvA -d 10000000 -t DP -f {} --output {} {}'.format(ref_faidx, vcf_mpileup_file, sort_bam_file)
    # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

    with unique_file_cm(ext='vcf') as vcf_mpileup_file:
        samtools_mpileup_with_defaults("--output", vcf_mpileup_file,
                                       "-f", ref_faidx,
                                       sort_bam_file
                                       )

    return vcf_mpileup_file
