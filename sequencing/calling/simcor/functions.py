import itertools

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