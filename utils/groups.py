__author__ = 'ofirr'

import math
from random import randint

def random_groups(groups, random_groups_keys, group_sizes_range=(10,24)):
    """create a dict representing groups with random sizes within range"""
    for group in random_groups_keys:
        groups[group]['size'] = randint(*group_sizes_range)
    return groups

def normalize_random_groups_to_sum(groups, random_groups_keys, random_groups_sum, group_sizes_range=(10,24)):
    """normalize sizes of random groups so that sum is fixed"""
    temp_sum = reduce(lambda x, y: x + y, [groups[group]['size'] for group in random_groups_keys])
    diff = random_groups_sum - temp_sum
    if diff == 0:
        return groups
    sign = int(math.copysign(1, diff))
    while abs(diff) > 0:
        group = random_groups_keys[randint(0, len(random_groups_keys)-1)]
        new_value = groups[group]['size'] + sign
        if new_value in range(*group_sizes_range):
            groups[group]['size'] += sign
            diff -= sign
    return groups

def random_groups_of_sum(groups, random_groups_keys, random_groups_sum, group_sizes_range=(10,24)):
    """get normalized dict of random sizes"""
    groups = random_groups(groups, random_groups_keys, group_sizes_range)
    return normalize_random_groups_to_sum(groups, random_groups_keys, random_groups_sum, group_sizes_range)