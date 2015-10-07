import re
from genomes.models import TargetType, Assembly, Target, Sequence, Chromosome
import hashlib

def add_target(name, target_type, assembly, chromosome, start_pos, end_pos, seq):
    """
    Checks input targets, validates against reference genome and saves target to DB
    returns (target, status_code, status_string)
    """
    try:
        target = Target.objects.get(start_pos=start_pos, end_pos=end_pos)
        return target, 0, "already exists"
    except Target.DoesNotExist:
        pass
    target = Target()
    target.name = name
    try:
        target.type = TargetType.objects.get(pk=target_type)
    except TargetType.DoesNotExist:
        return None, -1, 'Unknown target type: ' + target_type
    try:
        assembly_o = Assembly.objects.get(pk=assembly)
    except Assembly.DoesNotExist:
        return None, -1, 'Unknown assembly: ' + assembly
    try:
        target.chromosome = Chromosome.objects.get(name=chromosome, assembly=assembly_o)
    except Chromosome.DoesNotExist:
        return None, -1, 'Unknown chromosome: ' + chromosome
    target.start_pos = start_pos
    target.end_pos = end_pos
    if not re.match('^[ACTGactg]+$', seq.strip()):
        return None, -1, 'Corrupted sequence: ' + seq
    seq = seq.strip().upper()
    try:
        sequence = Sequence.objects.get(hash=hashlib.md5(seq).hexdigest())
    except Sequence.DoesNotExist:
        sequence = Sequence(length=len(seq), sequence=seq, hash=hashlib.md5(seq).hexdigest())
        sequence.save()
    target.referencevalue = sequence
    if target.chromosome.getdna(start_pos, end_pos) != sequence:
        try:
            target.start_pos, target.end_pos = target.chromosome.locate(target.start_pos, target.end_pos, sequence.sequence)
            try:
                target = Target.objects.get(start_pos=target.start_pos, end_pos=target.end_pos)
                return target, 0, "already exists"
            except Target.DoesNotExist:
                pass
        except ValueError:
            return None, -1, 'sequence does not match reference: %s\t%d\t%d' % (seq, start_pos, end_pos)
    target.save()
    return target, 1, "created"

# def getTabDelimitedValues(f):
#     """
#     Trivial tab-delimited parser
#     """
#     content = f.read()
#     f.close()
#     lines = []
#     for line in content.split('\n'):
#         lines.append(line.split('\t'))
#     return lines
#
# TARGETS_LINE_LENGTH = 9
#
#
# class TargetRow(object):
#     """
#     Define table columns
#     """
#     name = 0
#     type = 1
#     assembly = 2
#     chromosome = 3
#     start_pos = 4
#     end_pos = 5
#     fwd_primer = 6
#     rev_primer = 7
#     seq = 8
#
# def parseCommandsFromFile(file):
#     tabDelimitedValues = getTabDelimitedValues(file)
#     targets = []
#     for row in tabDelimitedValues:
#         if len(row) == TARGETS_LINE_LENGTH:
#             target = Target()
#             assert re.match('^[0-9]+$', row[TargetRow.type])
#             try:
#                 target.type = TargetType.objects.get(name = row[TargetRow.type])
#             except TargetType.DoesNotExist:
#                 print 'Unknown target type: ' + row[TargetRow.type]
#                 raise
#             assert re.match('^[0-9]+$', row[TargetRow.assembly])
#             try:
#                 target.assmbly = Assembly.objects.get(name = row[TargetRow.assembly])
#             except Assembly.DoesNotExist:
#                 print 'Unknown assembly: ' + row[TargetRow.assembly]
#                 raise
#             assert re.match('^chr[a-z]+$', row[TargetRow.chromosome])
#             target.chromosome = row[TargetRow.chromosome][3:]
#             # Check primers
#             assert re.match('^[0-9]+$', row[TargetRow.start_pos])
#             assert re.match('^[0-9]+$', row[TargetRow.end_pos])
#             assert re.match('^[ACTGactg]+$', row[TargetRow.fwd_primer])
#             assert re.match('^[ACTGactg]+$', row[TargetRow.rev_primer])
#             assert re.match('^[ACTGactg]+$', row[TargetRow.seq])
#             primer_fwd = row[TargetRow.fwd_primer]
#             assert row[TargetRow.seq][:len(primer_fwd)] == primer_fwd
#             primer_rev = complement(row[TargetRow.rev_primer]).reverse()
#             assert row[TargetRow.seq][-len(primer_rev):] == primer_rev
#             p = re.compile(r'[actg]+')
#             target_start = p.search(row[TargetRow.seq]).start()
#             target_end = p.search(row[TargetRow.seq]).end()
#
#             primer_forward = Primer()
#             primer_reverse = Primer()
#             #pf
#             primer_forward.start_pos = row[TargetRow.start_pos] - target_start
#             primer_forward.end_pos = row[TargetRow.start_pos] - target_start + len(row[TargetRow.fwd_primer])
#             try:
#                 Sequence.objects.get()
#             sequence = row[TargetRow.fwd_primer]
#             #pr
#             start_index = row[TargetRow.end_pos] - target_start
#             #target
#
#             try:
#                 primer_fwd = Primer.objects.get()
#
#