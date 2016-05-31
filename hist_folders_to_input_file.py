import glob
import os
import math
import sys
import argparse
import csv
import gzip
from frogress import bar as tqdm
import frogress
from collections import defaultdict
from itertools import combinations
from django.core.management import setup_environ
import sys
sys.path.append('/home/ofirr/CLineage/')
from clineage import settings
setup_environ(settings)
from linapp.models import Microsatellite, Panel, Target, TargetEnrichment, Chromosome, Plate
from scipy import stats
from order.hist import Histogram


def parse(path):
    cells = []
    with open(path) as f:
        rdr = csv.DictReader(f, dialect='excel-tab')
        for row in rdr:
            cell_row = {'path':row['histpath']}
            cells.append(cell_row)
    if not cells:
        print("""
no cells, file format should contain headers and formatted as follows:
histpath
/net/mraid11/export/data/dcsoft/home/LINEAGE/Miseq/MISEQR13/fastq/Output/A4-A2-cE2-AAR9-well-E5-BC-40_S37_L001_R1_001.hist
        """)
        raise
    return cells


def read_hist_file(path):
    hlines = open(path).readlines()[1:]
    return {line.split('\t',1)[0] : Histogram(list(map(int,line.split('\t')[1:-1])), normalize=False, nsamples=None, trunc=False, cut_peak=False) for line in hlines}


def read_hists(cells, verbose=True):
    return {cell_row['path'] : read_hist_file(cell_row['path']) for cell_row in tqdm(cells)}


def format_data(cells, loci_by_repeat_type):
    s = '\t'.join(['loc_name', 'cell_path'])
    s += '\t' + ' '.join([str(i+3) for i in range(148)])
    s += '\r\n'
    for cell in tqdm(cells):
        for repeat_type in ['AC']:
            for loci in loci_by_repeat_type[repeat_type].keys():
                if not loci_by_repeat_type[repeat_type][loci][cell['path']]:
                    continue
                if not loci_by_repeat_type[repeat_type][loci][cell['path']].nsamples:
                    continue
                row = '\t'.join([loci, cell['path']])
                row += '\t'
                for k in [(i+3) for i in range(148)]:
                    row += '{} '.format(loci_by_repeat_type[repeat_type][loci][cell['path']][k])
                s += row + '\r\n'
    return s


if '__main__' == __name__:
    parser = argparse.ArgumentParser(description='Generates histogram pairs')
    parser.add_argument('-i', '--input', type=str, dest='cells_table_file', help='path for tab delimited histpath-dup_id file')
    parser.add_argument('-o', '--output', type=str, dest='output', help='path for output file')
    parser.add_argument('-p', '--panel', type=str, dest='panel_name', default='aar9_panel', help='name of primers panel used for amplification')
    parser.add_argument('-a', '--assembly', type=str, dest='assembly_name', default='hg19', help='name of assembly (mm9/hg19...)')
    args = parser.parse_args()
    cells_table_file = args.cells_table_file
    output_file = args.output
    panel_name = args.panel_name
    assembly_name = args.assembly_name
    print('parsing input file')
    cells = parse(cells_table_file)
    print('parsing hists')
    cells_hists = read_hists(cells)
    loci = {}
    hists = []
    for cell, locs in cells_hists.items():
        for loc in locs.keys():
            loci.setdefault(loc, dict())[cell] = locs[loc]
            hists.append(locs[loc])
    print('retrieving contextual data from DB')
    X_chr = Chromosome.objects.get(assembly__friendly_name=assembly_name, name='X')
    try:
        panel = Panel.objects.get(name=panel_name)
    except Panel.DoesNotExist:
        print('No such panel: {}'.format(panel_name))
        exit()
    
    mss = Microsatellite.objects.filter(chromosome=X_chr).filter(primer_pair__in = panel.targets.all())
    # mss = Microsatellite.objects.exclude(chromosome=X_chr).filter(primer_pair__in = panel.targets.all())
    loci_by_repeat_type = defaultdict(lambda: defaultdict(dict))
    for ms in mss.filter(name__in=[key.split(":")[1] for key in loci.keys()]):
        for key in loci.keys():
            if ms.name == key.split(":")[1]:
                loci_by_repeat_type[ms.repeat_unit_type][key] = loci[key]
    print('formatting')
    s = format_data(cells, loci_by_repeat_type)
    print()
    print('writing file')
    f = gzip.open(output_file, 'wb')
    f.write(s)
    f.close()
    
