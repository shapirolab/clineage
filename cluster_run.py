import os
os.chdir('/home/ofirr/s/Ofir/hcc/hist_calling')
from order.utils.output_formatters import load_or_create_calling
from order.utils.output_formatters import generate_output_file, load_or_create_simulations_file, save_calling_file
import numpy
from frogress import bar
from IPython.parallel import Client
from itertools import izip, repeat
import itertools

normalize = True
truncate = False
cutpeak = False
trim_extremes = False
# output_file = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Calling/NSR3_AC_nonX_bin_1a.tab.gz'
output_file = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Miseq/ML76/Unified/fastq/Output/Calling/ML76_AC_X_mat_1a.tab.gz'
# input_file = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Calling/input_hists_AC_nonX.tab.gz'
input_file = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Miseq/ML76/Unified/fastq/Output/Calling/input_hists_AC_X.tab.gz'
reads_threshold = 30
score_threshold = 999.9
cycles_threshold = 60
min_cycles = 20
shift_margins = 8
max_alleles = 1
sample_depth = 10000
max_ms_length = 50
max_distance_from_median = 3
meth = 'con'
verbose = True
sim_method = 'mat'
SIMULATED_HISTS_PATH = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Miseq/ML76/Unified/fastq/Output/Calling/ac_mat_1a_sim_hists_15o.pickle'
SIGNALS_CALLING_PATH = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Miseq/ML76/Unified/fastq/Output/Calling/ML76_AC_X_mat_1a.pickle'

# AC #
ups = [numpy.poly1d([0.00026892, -0.00205025])]
dws = [numpy.poly1d([0.00191615, -0.01174076]),
       numpy.poly1d([0.00027444, -0.00220836]),
       numpy.poly1d([0.0001768, -0.00199328])]

sim_hists = load_or_create_simulations_file(SIMULATED_HISTS_PATH,
                                            method=sim_method,
                                            max_cycles=cycles_threshold,
                                            ups=ups,
                                            dws=dws,
                                            max_ms_length=max_ms_length,
                                            sample_depth=sample_depth,
                                            max_alleles=max_alleles,
                                            normalize=normalize,
                                            truncate=truncate,
                                            cut_peak=cutpeak,
                                            trim_extremes=trim_extremes)

calling = load_or_create_calling(SIGNALS_CALLING_PATH)

rc = Client()
dview = rc[:]
with dview.sync_imports():
    import os
dview.execute("os.chdir('/home/ofirr/s/Ofir/hcc/hist_calling/')")
with dview.sync_imports():
    from order.calling import call_multi_hist, inflate_index, flatten_index, uncalled_inputs

lview = rc.load_balanced_view()
lview.block = False
cpr = dview.use_cloudpickle()


def grouper(n, iterable):
    it = iter(iterable)
    while True:
       chunk = tuple(itertools.islice(it, n))
       if not chunk:
           return
       yield chunk


def helper(input_tup):
    results = []
    uncalled_inputs_group, extras = input_tup
    for loc, cell, row_hist in uncalled_inputs_group:
        flat_sim_hists, kwargs = extras
        sim_hists = inflate_index(flat_sim_hists)
        res = call_multi_hist(row_hist, sim_hists, **kwargs)
        results.append((loc, cell, row_hist, res))
    return results

kwargs = {'method': meth,
          'score_threshold': score_threshold,
          'min_cycles': min_cycles,
          'max_cycles': cycles_threshold,
          'nsamples': None,
          'normalize': normalize,
          'truncate': truncate,
          'cut_peak': cutpeak,
          'trim_extremes': trim_extremes,
          'shift_margins': shift_margins,
          'max_distance_from_median': max_distance_from_median,
          'max_alleles': max_alleles,
          'max_ms_length': max_ms_length}


uncalled_inputs_groups = grouper(10000, uncalled_inputs(input_file, calling, reads_threshold=reads_threshold))
flat_sim_hists = flatten_index(sim_hists)
inputs_generator = izip(uncalled_inputs_groups, repeat((flat_sim_hists, kwargs)))

for results in bar(lview.map(helper, inputs_generator)):
    for result in results:
        loc, cell, row_hist, res = result
        calling[loc][cell] = res


save_calling_file(calling, SIGNALS_CALLING_PATH)
generate_output_file(input_file,
                     output_file,
                     calling,
                     reads_threshold=reads_threshold,
                     score_threshold=score_threshold,
                     verbose=verbose)
