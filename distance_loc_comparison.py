__author__ = 'ofirr'
import os
os.chdir('/home/ofirr/s/Ofir/hcc/hist_calling')
from order.utils.output_formatters import generate_output_file, load_or_create_simulations_file, save_calling_file
import numpy
from frogress import bar
from itertools import izip, repeat
import itertools
from order.utils.parsers import parse_input_file

normalize = True
truncate = False
cutpeak = False
trim_extremes = False
output_file = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Calling/NSR3_AC_nonX_mat_2a_prop.tab.gz'
input_file = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Calling/input_hists_AC_nonX.tab.gz'
reads_threshold = 30
score_threshold = 999.9
cycles_threshold = 50
min_cycles = 30
shift_margins = 16
max_alleles = 2
sample_depth = 10000
min_ms_length = 5
max_ms_length = 35
max_distance_from_median = 8
# length_sensitivity = 0.08
# diff_sensetivity = 0.5
length_sensitivity = 0.4
diff_sensetivity = 1.5
# length_sensitivity = 0.21
# diff_sensetivity = 0.65
# length_sensitivity = 0.8
# diff_sensetivity = 1.3
# length_sensitivity = 0.01
# diff_sensetivity = 0.0
meth = 'con'
verbose = True
sim_method = 'mat'
SIMULATED_HISTS_PATH = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Calling/ac_mat_2a_prop_100_sim_hists_with_replacement.pickle'
SIGNALS_CALLING_PATH = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Calling/NSR3_AC_nonX_mat_2a_prop.pickle'


class hashable_poly1d(numpy.poly1d):
    def __hash__(self):
        return hash((tuple(self.coeffs), self.order, self.variable))


ups = [hashable_poly1d([0.00026892, -0.00205025])]
dws = [hashable_poly1d([0.00191615, -0.01174076]),
       hashable_poly1d([0.00027444, -0.00220836]),
       hashable_poly1d([0.0001768, -0.00199328])]

sim_hists = load_or_create_simulations_file(SIMULATED_HISTS_PATH,
                                            method=sim_method,
                                            min_cycles=min_cycles,
                                            max_cycles=cycles_threshold,
                                            ups=ups,
                                            dws=dws,
                                            min_ms_length=min_ms_length,
                                            max_ms_length=max_ms_length,
                                            sample_depth=sample_depth,
                                            max_alleles=max_alleles,
                                            normalize=normalize,
                                            truncate=truncate,
                                            cut_peak=cutpeak,
                                            trim_extremes=trim_extremes)

# for seeds in sim_hists:
#     for seeds_and_proportions in sim_hists[seeds]:
#         for cycle in sim_hists[seeds][seeds_and_proportions]:
#             sim_hists[seeds][seeds_and_proportions][cycle].normalize()


from collections import defaultdict
calling = defaultdict(lambda: defaultdict(dict))
from order.calling import call_multi_hist


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
        # flat_sim_hists, kwargs = extras
        sim_hists, kwargs = extras
        # sim_hists = inflate_index(flat_sim_hists)
        res = call_multi_hist(row_hist, sim_hists, proportional=True, **kwargs)
        #res = None
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
          'max_ms_length': max_ms_length,
          'length_sensitivity': length_sensitivity,
          'diff_sensetivity': diff_sensetivity}

H1_map = {
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC612-cellID2535-Tamir_S54_L001_R1_001.hist': (2535, 'A7_A5_cA1'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC636-cellID2536-Tamir_S86_L001_R1_001.hist': (2536, 'A7_A5_cB1'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC633-cellID2537-Tamir_S77_L001_R1_001.hist': (2537, 'A7_A5_cC1'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC632-cellID2539-Tamir_S59_L001_R1_001.hist': (2539, 'A7_A5_cE1'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC630-cellID2540-Tamir_S53_L001_R1_001.hist': (2540, 'A7_A5_cF1'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC652-cellID2541-Tamir_S83_L001_R1_001.hist': (2541, 'A7_A5_cG1'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC658-cellID2542-Tamir_S62_L001_R1_001.hist': (2542, 'A7_A5_cA2'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC592-cellID2543-Tamir_S19_L001_R1_001.hist': (2543, 'A7_A5_cB2'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC626-cellID2544-Tamir_S48_L001_R1_001.hist': (2544, 'A7_A5_cD2'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC614-cellID2545-Tamir_S58_L001_R1_001.hist': (2545, 'A7_A5_cE2'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC669-cellID2546-Tamir_S49_L001_R1_001.hist': (2546, 'A7_A5_cF2'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC627-cellID2547-Tamir_S56_L001_R1_001.hist': (2547, 'A7_A5_cB3'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC648-cellID2548-Tamir_S70_L001_R1_001.hist': (2548, 'A7_A5_cD3'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC670-cellID2549-Tamir_S41_L001_R1_001.hist': (2549, 'A7_A5_cE3'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC651-cellID2550-Tamir_S10_L001_R1_001.hist': (2550, 'A7_A5_cH3'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC602-cellID2551-Tamir_S31_L001_R1_001.hist': (2551, 'A6_B2_cB4'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC650-cellID2552-Tamir_S68_L001_R1_001.hist': (2552, 'A6_B2_cC4'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC606-cellID2553-Tamir_S38_L001_R1_001.hist': (2553, 'A6_B2_cD4'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC668-cellID2554-Tamir_S73_L001_R1_001.hist': (2554, 'A6_B2_cE4'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC605-cellID2555-Tamir_S35_L001_R1_001.hist': (2555, 'A6_B2_cF4'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC601-cellID2556-Tamir_S30_L001_R1_001.hist': (2556, 'A6_B2_cG4'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC619-cellID2557-Tamir_S76_L001_R1_001.hist': (2557, 'A6_B2_cH4'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC645-cellID2558-Tamir_S39_L001_R1_001.hist': (2558, 'A6_B2_cA5'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC625-cellID2559-Tamir_S51_L001_R1_001.hist': (2559, 'A6_B2_cB5'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC635-cellID2561-Tamir_S84_L001_R1_001.hist': (2561, 'A6_B2_cD5'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC657-cellID2562-Tamir_S17_L001_R1_001.hist': (2562, 'A6_B2_cE5'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC654-cellID2563-Tamir_S60_L001_R1_001.hist': (2563, 'A6_B2_cF5'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC600-cellID2564-Tamir_S29_L001_R1_001.hist': (2564, 'A6_B2_cG5'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC637-cellID2565-Tamir_S57_L001_R1_001.hist': (2565, 'A6_B2_cH5'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC591-cellID2566-Tamir_S18_L001_R1_001.hist': (2566, 'A6_B2_cA6'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC629-cellID2566-Tamir_S42_L001_R1_001.hist': (2566, 'A6_B2_cA6'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC641-cellID2567-Tamir_S74_L001_R1_001.hist': (2567, 'A6_B2_cB6'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC584-cellID2568-Tamir_S5_L001_R1_001.hist': (2568, 'A6_B2_cD6'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC639-cellID2569-Tamir_S78_L001_R1_001.hist': (2569, 'A6_B2_cE6'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC604-cellID2570-Tamir_S34_L001_R1_001.hist': (2570, 'A6_B2_cF6'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC608-cellID2571-Tamir_S45_L001_R1_001.hist': (2571, 'A6_B2_cH6'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC620-cellID2572-Tamir_S81_L001_R1_001.hist': (2572, 'A6_D1_cA7'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC583-cellID2573-Tamir_S4_L001_R1_001.hist': (2573, 'A6_D1_cC7'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC667-cellID2574-Tamir_S46_L001_R1_001.hist': (2574, 'A6_D1_cD7'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC587-cellID2575-Tamir_S9_L001_R1_001.hist': (2575, 'A6_D1_cE7'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC593-cellID2575-Tamir_S22_L001_R1_001.hist': (2575, 'A6_D1_cE7'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC596-cellID2577-Tamir_S25_L001_R1_001.hist': (2577, 'A6_D1_cG7'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC589-cellID2578-Tamir_S14_L001_R1_001.hist': (2578, 'A6_D1_cH7'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC590-cellID2579-Tamir_S16_L001_R1_001.hist': (2579, 'A6_D1_cA8'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC617-cellID2580-Tamir_S67_L001_R1_001.hist': (2580, 'A6_D1_cB8'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC578-cellID2581-Tamir_S1_L001_R1_001.hist': (2581, 'A6_D1_cD8'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC585-cellID2582-Tamir_S6_L001_R1_001.hist': (2582, 'A6_D1_cE8'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC631-cellID2583-Tamir_S69_L001_R1_001.hist': (2583, 'A6_D1_cF8'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC597-cellID2584-Tamir_S26_L001_R1_001.hist': (2584, 'A6_D1_cG8'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC664-cellID2585-Tamir_S43_L001_R1_001.hist': (2585, 'A6_D1_cH8'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC613-cellID2586-Tamir_S55_L001_R1_001.hist': (2586, 'A7_H9_cB9'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC665-cellID2587-Tamir_S79_L001_R1_001.hist': (2587, 'A7_H9_cD9'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC634-cellID2588-Tamir_S66_L001_R1_001.hist': (2588, 'A7_H9_cE9'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC649-cellID2589-Tamir_S44_L001_R1_001.hist': (2589, 'A7_H9_cF9'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC616-cellID2591-Tamir_S65_L001_R1_001.hist': (2591, 'A7_A5_F1_cB1'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC661-cellID2595-Tamir_S36_L001_R1_001.hist': (2595, 'A7_A5_F1_cG1'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC655-cellID2599-Tamir_S63_L001_R1_001.hist': (2599, 'A7_A5_F1_cF2'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC595-cellID2604-Tamir_S24_L001_R1_001.hist': (2604, 'A7_A5_F1_cC3'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC666-cellID2606-Tamir_S11_L001_R1_001.hist': (2606, 'A7_A5_F1_cE3'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC581-cellID2609-Tamir_S2_L001_R1_001.hist': (2609, 'A7_A5_F1_cH3'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC611-cellID2610-Tamir_S52_L001_R1_001.hist': (2610, 'A6_D1_E7_cB4'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC662-cellID2611-Tamir_S85_L001_R1_001.hist': (2611, 'A6_D1_E7_cC4'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC599-cellID2612-Tamir_S27_L001_R1_001.hist': (2612, 'A6_D1_E7_cD4'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC588-cellID2613-Tamir_S12_L001_R1_001.hist': (2613, 'A6_D1_E7_cE4'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC586-cellID2614-Tamir_S8_L001_R1_001.hist': (2614, 'A6_D1_E7_cF4'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC607-cellID2615-Tamir_S40_L001_R1_001.hist': (2615, 'A6_D1_E7_cG4'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC622-cellID2616-Tamir_S89_L001_R1_001.hist': (2616, 'A6_D1_E7_cH4'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC582-cellID2617-Tamir_S3_L001_R1_001.hist': (2617, 'A6_D1_E7_cA5'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC621-cellID2618-Tamir_S88_L001_R1_001.hist': (2618, 'A6_D1_E7_cB5'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC646-cellID2619-Tamir_S33_L001_R1_001.hist': (2619, 'A6_D1_E7_cC5'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC644-cellID2620-Tamir_S75_L001_R1_001.hist': (2620, 'A6_D1_E7_cE5'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC618-cellID2621-Tamir_S72_L001_R1_001.hist': (2621, 'A6_D1_E7_cF5'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC653-cellID2622-Tamir_S71_L001_R1_001.hist': (2622, 'A6_D1_E7_cH5'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC628-cellID2623-Tamir_S64_L001_R1_001.hist': (2623, 'A6_D1_E7_cA6'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC743-cellID2624-Tamir_S87_L001_R1_001.hist': (2624, 'A6_D1_E7_cB6'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC711-cellID2625-Tamir_S28_L001_R1_001.hist': (2625, 'A6_D1_E7_cC6'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC694-cellID2626-Tamir_S13_L001_R1_001.hist': (2626, 'A6_D1_E7_cD6'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC678-cellID2627-Tamir_S37_L001_R1_001.hist': (2627, 'A6_D1_E7_cE6'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC699-cellID2627-Tamir_S21_L001_R1_001.hist': (2627, 'A6_D1_E7_cE6'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC697-cellID2628-Tamir_S32_L001_R1_001.hist': (2628, 'A6_D1_E7_cF6'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC759-cellID2630-Tamir_S80_L001_R1_001.hist': (2630, 'A6_D1_E7_cH6'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC713-cellID2631-Tamir_S23_L001_R1_001.hist': (2631, 'A6_D1_E7_cA7'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC705-cellID2632-Tamir_S91_L001_R1_001.hist': (2632, 'A6_D1_E7_cB7'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC690-cellID2633-Tamir_S20_L001_R1_001.hist': (2633, 'A6_D1_E7_cC7'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC734-cellID2633-Tamir_S82_L001_R1_001.hist': (2633, 'A6_D1_E7_cC7'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC677-cellID2635-Tamir_S90_L001_R1_001.hist': (2635, 'A6_D1_E7_cE7'),
'/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Output/BC757-cellID2635-Tamir_S47_L001_R1_001.hist': (2635, 'A6_D1_E7_cE7'),
}

def uncalled_inputs_for_loc(input_file, calling, H1_map, locs_filter, reads_threshold=0):
    for loc, cell, row_hist in parse_input_file(input_file):
        if loc not in locs_filter:
            continue
        if calling[loc][cell]:
            continue
        if cell not in H1_map:
            continue
        if sum(row_hist.values()) < reads_threshold:
            continue
        yield loc, cell, row_hist

#loc = '29892:5_65361420_65361487'  # 23, 24
#loc = '29832:16_51872847_51872908'  # 13, 15
#loc = '37926:7_76889257_76889314'  #21, 23
#loc = '29927:9_85583222_85583285'  # 25, 28
#loc = '59820:1_87083880_87083935'  # 22, 27
#loc = '29937:12_11902223_11902284'  # 17, 27
#loc = '29910:7_92001265_92001326'  # 27,28


# for meth in ['con', 'kso', 'ksa', 'cor', 'sub', 'sp', 'chi', 'pr', 'emd', 'aks', 'ks', 'a2ks']:
for meth in ['con']:
    for loc in ['29927:9_85583222_85583285', '29832:16_51872847_51872908', '37926:7_76889257_76889314', '59820:1_87083880_87083935', '29937:12_11902223_11902284', '29892:5_65361420_65361487']:
    # for loc in ['29892:5_65361420_65361487']:
        kwargs['method'] = meth
        calling = defaultdict(lambda: defaultdict(dict))
        uncalled_inputs_groups = list(grouper(1, uncalled_inputs_for_loc(input_file, calling, H1_map, locs_filter=[loc], reads_threshold=reads_threshold)))
        inputs_generator = list(izip(uncalled_inputs_groups, repeat((sim_hists, kwargs))))
        for input_params in inputs_generator:
            results = helper(input_params)
            for result in results:
                loc, cell, row_hist, res = result
                calling[loc][cell] = res
                print meth, calling[loc][cell]['seeds_and_proportions'], calling[loc][cell]['score'], calling[loc][cell]['reads'], loc, cell
                #print meth, loc, cell