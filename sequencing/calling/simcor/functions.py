import itertools
import decimal

def get_far_apart_highest_peaks(hist, k, d):
    #hist is the histogram
    #k is the allele_number
    #d is the minimum distance between the allele
    hs = sorted(hist._hist.items(), key=lambda hkey: hkey[1], reverse=True)
    for allele_number in range(k):
        if allele_number >= len(hs):
            break
        for g in range(1, d):
            if [item for item in hs if item[0] == hs[allele_number][0]+g]:
                hs.remove([item for item in hs if item[0] == hs[allele_number][0]+g][0])
            if [item for item in hs if item[0] == hs[allele_number][0]-g]:
                hs.remove([item for item in hs if item[0] == hs[allele_number][0]-g][0])
    seeds = [x for x, y in hs[:k]]
    return seeds


def seeds_search_range( peaks, max_distance_between_peaks, max_ms_length ):
    search_range= itertools.product(
            *[
                range(
                    max(1, peak-max_distance_between_peaks),
                    min(max_ms_length, peak+max_distance_between_peaks+1)
                ) for peak in peaks
            ])
    yield from search_range



"""
line of events:
 get_guide
 separate_alleles
 done
"""

def get_sorted_seeds(seeds_and_proportions):
    gs = zip (*seeds_and_proportions)
    return sorted(list(gs)[0])


def get_proportions(seeds_and_proportions):
    gs = zip (*seeds_and_proportions)
    return sorted(list(gs)[1])

def get_guide(loc_seeds):
    guide=[]
    guide_points=[]
    for cell in loc_seeds():
        cell_seeds=loc_seeds[cell]['seeds_and_proportions']
        number_of_alleles = len(cell_seeds)
        if number_of_alleles > 1:
            wanted_proportion = 1 / number_of_alleles
            if guide:
                if( abs(decimal.Decimal(get_proportions(guide)[0]) - decimal.Decimal(wanted_proportion )) >
                    abs (decimal.Decimal(get_proportions(cell_seeds)[0]) - decimal.Decimal(wanted_proportion ))
                  ) :
                        guide=cell_seeds
            if not guide:
                guide=cell_seeds
    return guide

def extract_seeds(loc_seeds):
    extracted_seeds={}
    for cell in loc_seeds:
        cell_seeds=loc_seeds[cell]['seeds_and_proportions']
        if len(cell_seeds)>1:
            extracted_seeds[cell] = get_sorted_seeds(cell_seeds)
        else:
            extracted_seeds[cell] = [seed for seed in cell_seeds]
    return extracted_seeds

def separate_alleles(loc_seeds, guide_points, loc_name='loci'):
    extracted_seed= extract_seeds(loc_seeds)
    separated_loc={}
    guide_points_allels=[]
    for allele in guide_points:
        guide_points_allels.append('{}_{}'.format(loc_name,chr(guide_points.index(allele)+ord('a'))))

    for cell in extracted_seeds:
        for seed in extracted_seeds[cell]:
            closest_to_guide_index=min(guide_points, key=lambda x:abs(x-seed))
            if (abs(seed-closest_to_guide_index)<=3):
                separated_loc.setdefault(guide_points_allels[guide_points.index(closest_to_guide_index)],dict())[cell]=seed
    return separated_loc




