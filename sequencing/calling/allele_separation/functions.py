import decimal

"""
line of events:
 get_guide
 separate_alleles
 done
"""

def get_sorted_seeds(seeds_and_proportions):
    gs = zip(*seeds_and_proportions)
    return sorted(list(gs)[0])


def get_proportions(seeds_and_proportions):
    gs = zip(*seeds_and_proportions)
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
    extracted_seeds= extract_seeds(loc_seeds)
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
