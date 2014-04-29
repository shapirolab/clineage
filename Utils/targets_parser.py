
import re
from LinApp.models import TargetType, Assembly, Target, Primer, Sequence

def getTabDelimitedValues(f):
    """
    Trivial tab-delimited parser
    """
    content = f.read()
    f.close()
    lines = []
    for line in content.split('\n'):
        lines.append(line.split('\t'))
    return lines

TARGETS_LINE_LENGTH = 9


class TargetRow(object):
    """
    Define table columns
    """
    name = 0
    type = 1
    assembly = 2
    chromosome = 3
    start_pos = 4
    end_pos = 5
    fwd_primer = 6
    rev_primer = 7
    seq = 8

def parseCommandsFromFile(file):
    tabDelimitedValues = getTabDelimitedValues(file)
    targets = []
    for row in tabDelimitedValues:
        if len(row) == COMMANDS_LINE_LENGTH:
            target = Target()
            assert re.match('^[0-9]+$', row[TargetRow.type])
            try:
                target.type = TargetType.objects.get(name = row[TargetRow.type])
            except ObjectDoesNotExist:
                print 'Unknown target type: ' + row[TargetRow.type]
                raise
            assert re.match('^[0-9]+$', row[TargetRow.assembly])
            try:
                target.assmbly = Assembly.objects.get(name = row[TargetRow.assembly])
            except ObjectDoesNotExist:
                print 'Unknown assembly: ' + row[TargetRow.assembly]
                raise
            assert re.match('^chr[a-z]+$', row[TargetRow.chromosome])
            target.chromosome = row[TargetRow.chromosome][3:]
            # Check primers
            assert re.match('^[0-9]+$', row[TargetRow.start_pos])
            assert re.match('^[0-9]+$', row[TargetRow.end_pos])
            assert re.match('^[ACTGactg]+$', row[TargetRow.fwd_primer])
            assert re.match('^[ACTGactg]+$', row[TargetRow.rev_primer])
            assert re.match('^[ACTGactg]+$', row[TargetRow.seq])
            primer_fwd = row[TargetRow.fwd_primer])
            assert row[TargetRow.seq][:len(primer_fwd)] == primer_fwd
            primer_rev = complement(row[TargetRow.rev_primer]).reverse()
            assert row[TargetRow.seq][-len(primer_rev):] == primer_rev
            p = re.compile(r'[actg]+')
            target_start = p.search(row[TargetRow.seq]).start()
            target_end = p.search(row[TargetRow.seq]).end()
            
            primer_forward = Primer()
            primer_reverse = Primer()
            #pf
            primer_forward.start_pos = row[TargetRow.start_pos] - target_start
            primer_forward.end_pos = row[TargetRow.start_pos] - target_start + len(row[TargetRow.fwd_primer])
            try:
                Sequence.objects.get()
            sequence = row[TargetRow.fwd_primer]
            #pr
            start_index = row[TargetRow.end_pos] - target_start
            #target
            
            try:
                primer_fwd = Primer.objects.get()
            
            