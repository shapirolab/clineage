from lib_prep.models import Sequencing, Machine
import os
import re
import gzip

def

def parse_sample_sheet(path):
    try:
        csv_file = open(path,'rb')
        s = csv_file.read()
    except


def import_sequencing(seq, match=None):
    if match is not None:
        match_re = re.compile(match)
    ouput_path = '\\\\%s\\%s' % (seq.machine.ip, 'd\\Illumina\\MiSeqOutput')
    for f in sorted(os.listdir(ouput_path)):
        sequencing_directory = os.path.join(ouput_path, f)
        if os.path.isdir(sequencing_directory) and (match is None or match_re.search(f)):
            sample_sheet_path = '%s\\%s\\%s' % (f, 'Data\Intensities\BaseCalls','SampleSheet.csv')
            if os.path.isfile(sample_sheet_path):
                cells = parse_sample_sheet(sample_sheet_path)
            for cell_fastq_file in os.listdir('%s\\%s' % (f, 'Data\Intensities\BaseCalls')):
                if