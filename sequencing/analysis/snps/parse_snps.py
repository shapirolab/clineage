# __author__ = 'veronika'

from collections import defaultdict
import pickle
from misc.utils import unique_file_cm
from targeted_enrichment.planning.models import Target, SNP, Microsatellite
from targeted_enrichment.amplicons.models import Amplicon, UMITargetedAmplicon


base_map = {'A': 1,
            'C': 2,
            'G': 3,
            'T': 4,
            'a': 1,
            'c': 2,
            'g': 3,
            't': 4}


def print_to_file(file_name, rows):
    with open(file_name, 'w') as f:
        for row in rows:
            f.write(row)
    pass


def valid_call(vcf_loc, min_cover=5):
    if not vcf_loc:
        return False
    if vcf_loc['DP'] < min_cover:
        return False
    for i in vcf_loc['stats']:
        if '-nan' == i:
            return False
    return True


def get_locs_cells(snp_dict):
    locs = set()
    cells = set()
    for loc in snp_dict:
        for cell in snp_dict[loc]:
            locs.add(loc)
            cells.add(cell)
    return list(locs), list(cells)


def get_method(method):
    if method == 'read_count':
        return read_count  # TODO: change to num_reads_histogram
    if method == 'explicit_snps':
        return explicit_snp_calling
    if method == 'allele_count':
        return allele_count
    if method == 'full':
        return snp_calling
    if method == 'snp_defind':
        return snp_defind
    raise NotImplemented


def read_count(vcf_loc, score_threshold, min_cover=5):
    if not valid_call(vcf_loc, min_cover):
        return 'NaN'
    return vcf_loc['DP']


def snp_defind(vcf_loc, score_threshold, min_cover=5):
    if not valid_call(vcf_loc, min_cover):
        return 'NaN'
    call = snp_calling(vcf_loc, score_threshold, min_cover)
    soreted_bases = sorted(call.items(), key=lambda x: x[1], reverse=True)
    base, prob = soreted_bases[0]

    return base


def explicit_snp_calling(vcf_loc, score_threshold, min_cover=5):
    if not valid_call(vcf_loc, min_cover):
        return 'NaN'
    call = snp_calling(vcf_loc, score_threshold, min_cover)
    soreted_bases = sorted(call.items(), key=lambda x: x[1], reverse=True)
    alleles = []
    for base, prob in soreted_bases:
        # mono allelic case, 1.0 would be the most strict, all reads cluster around that single allele
        if prob >= score_threshold or (1.5 - score_threshold >= prob >= score_threshold - 0.4999):
            if vcf_loc['SNP_defined']:
                alleles.append(base)
    return alleles


def allele_count(vcf_loc, score_threshold, min_cover=5):
    if not valid_call(vcf_loc, min_cover):
        return 'NaN'
    return len((vcf_loc, score_threshold, min_cover))


def snp_calling(vcf_loc, score_threshold, min_cover=5):
    if not valid_call(vcf_loc, min_cover):
        return 'NaN'
    call = dict()
    stats = [float(i) for i in vcf_loc['stats']]
    call[vcf_loc['base']] = stats[0]
    for i, st in enumerate(stats[1:]):
        if st > 0 and len(vcf_loc['modification']) >= i + 1:
            if str(type(vcf_loc['modification'])) == "<class 'Bio.SeqRecord.SeqRecord'>":
                call[vcf_loc['modification'].seq] = st
            else:
                call[vcf_loc['modification'][i].sequence] = st
    return call


def convert_snp_dict_to_table(snp_dict, score_threshold, method, min_cover=5):
    rows = []
    locs, cells = get_locs_cells(snp_dict)
    calling_method = get_method(method)
    rows.append('Loci\t' + '\t'.join(cells) + '\n')
    for loc in locs:
        if method == 'explicit_snps':
            if snp_dict[loc][list(snp_dict[loc].keys())[0]]['SNP_defined']:
                row = '{}.{}\t'.format(*loc)
            else:
                continue
        else:
            if snp_dict[loc][list(snp_dict[loc].keys())[0]]['SNP_defined']:
                row = '{}.{}*\t'.format(*loc)
            else:
                row = '{}.{}\t'.format(*loc)

        for cell in cells:
            if cell in snp_dict[loc]:
                call_num = calling_method(snp_dict[loc][cell], score_threshold, min_cover)
                row += '{}\t'.format(str(call_num))
            else:
                row += 'NaN\t'
        row += '\n'
        rows.append(row)
    return rows


###############################################################
def get_locs_amplicons(locs):
    locs_amplicons = {}
    amplicons_locs = {list()}
    for loc in locs:
        amplicon = loc[0]
        locs_amplicons[loc] = amplicon
        amplicons_locs[amplicon].append(loc)
    return amplicons_locs, locs_amplicons


def get_amplicons_reads_map(snp_dict):
    locs, cells = get_locs_cells(snp_dict)
    amplicons_locs, locs_amplicons = get_locs_amplicons(locs)
    amplicons_reads = defaultdict(lambda: defaultdict(int))
    for cell in cells:
        for amplicon in amplicons_locs:
            reads = 0
            for loc in amplicons_locs[amplicon]:
                if not valid_call(snp_dict[loc][cell], min_cover=0):
                    continue
                reads += snp_dict[loc][cell]['DP']
            call_num = reads/len(amplicons_locs[amplicon])
            amplicons_reads[amplicon][cell] = call_num
    return amplicons_reads


def snp_is_ms(snp):
    amp = Amplicon.objects.get(id=snp.CHROM)
    slices = amp.slice.contains.all()
    for sl in slices:
        if not sl.target_set.count():
            continue
        for tr in Target.objects.filter(slice=sl):
            try:
                tr.microsatellite
                return True
            except Microsatellite.DoesNotExist:
                pass
    return False


def get_rel_pos(amplicon, snp):
    """
    Returns the relative position of a SNP on an amplicon.
              |--------------------------------------------------------|
     amplicon.slice: |------------------------------------------|
           UMI|PRIMER--------------T-----------------------------PRIMER|UMI
    """
    if amplicon.slice.start_pos == amplicon.left_ugs.slice.start_pos:
        rel_pos = amplicon.slice.relative_pos(snp.slice.chromosome, snp.slice.start_pos) + amplicon.umi_length + 1
    else:
        rel_pos = amplicon.slice.relative_pos(snp.slice.chromosome, snp.slice.start_pos) + len(amplicon.left_margin)
    return rel_pos


def get_amplicon_snps(amplicon):
    for sl in amplicon.slice.contains.all():
        for tr in Target.objects.filter(slice=sl):
            try:
                tr.snp
                yield tr
            except SNP.DoesNotExist:
                pass


def retrieve_explicit_snps_positions(vcf_snp):
    amplicon = Amplicon.objects.get(id=vcf_snp.CHROM).subclass
    for snp in get_amplicon_snps(amplicon):
        yield get_rel_pos(amplicon, snp)


def create_snp_file(snps_dict):
    with unique_file_cm(ext='pickle') as handle:
        pickle.dump(snps_dict, open(handle, 'wb'))

    return handle
