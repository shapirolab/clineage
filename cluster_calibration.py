from cloud.serialization.cloudpickle import dumps, loads
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from collections import Counter
from scipy import optimize
import os
# %matplotlib inline

os.chdir('/home/ofirr/s/Ofir/hcc/hist_calling/')
from order.hist import Histogram
from order.preprocessing import generate_hist
from order.calibration.score import optimize_across_lengths
from order.calibration.score import param_list_to_polynomes


def plot_result(ret, protocols, syn_hists, repeat_unit):
    params = ret[5]['x']
    score = ret[5]['fun']
    ups, dws = param_list_to_polynomes(params)
    base = int(len(syn_hists[repeat_unit])**0.5)+1
    # plt.close('all')
    # fig = plt.figure(figsize=(10,10))
    # gs = gridspec.GridSpec(base, base)
    for repeat_number in syn_hists[repeat_unit]:
        for c, t in protocols:
            model_hist = generate_hist(repeat_number, c, sim, ups=ups, dws=dws) + repeat_number
    #         print repeat_number, c, sim, ups, dws
    #         print model_hist
            real_hist = syn_hists[repeat_unit][repeat_number][t][0]
            model_hist.sq_normalize()
            real_hist.sq_normalize()
            # fig = plt.figure(figsize=(10,7))
            plt.xlim(0,repeat_number+20)
            x = np.arange(0,repeat_number+100)
            plt.scatter(x, [real_hist._hist[i] for i in x], s=10)
            plt.plot(x, [real_hist._hist[i] for i in x], label='real {} {} {} {}'.format(repeat_number, repeat_unit, c, real_hist.nsamples))
            plt.scatter(x, [model_hist._hist[i] for i in x], s=10)
            plt.plot(x, [model_hist._hist[i] for i in x], label='model {} {} {}'.format(repeat_number, repeat_unit, t))
            plt.legend()
            plt.show()


with open('experimental_data/synthetic_hists.pickle', 'rb') as f:
    syn_hists = loads(f.read())

c1, c2 = 14, 45
repeat_unit = 'AC'
t1 = 'PCR2only'
t2 = 'PCR1-2'
hist_pairs = [(repeat_number, syn_hists[repeat_unit][repeat_number][t1][0]-repeat_number, syn_hists[repeat_unit][repeat_number][t2][0]-repeat_number) for repeat_number in syn_hists[repeat_unit]]
#hist_pairs = get_transposed_hists_pairs(filename='experimental_data/hist_by_ms_len_as_0_sum.tab', lengths=[30])
alg = 'con'
# optimizer_method = "L-BFGS-B"
optimizer_method = "SLSQP"
bounds = [(0.000, 0.0001), (-0.01, 0.01), (-0.1, 0.1),
          (0.000, 0.0001), (-0.01, 0.01), (-0.1, 0.1),]
optimizer_options = {'eps': 1e-5, 'disp': True}
initial_guess = [0, 0, 0,
                 0, 0, 0,]
iterations = 100
sim = 'bon'
ret = optimize_across_lengths((alg, sim, optimizer_method, hist_pairs, (c1, c2), bounds, initial_guess, iterations, optimizer_options))
