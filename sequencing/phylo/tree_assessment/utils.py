import random
import numpy
from itertools import combinations


def subdict(bigdict, subset_keys):
    intersection_keys = bigdict.keys() & set(subset_keys)
    return {k: bigdict[k] for k in intersection_keys}


def random_choose(l, k):
    "Randomized itertools.combinations(iterable, r)"
    vs = list(l)
    random.shuffle(vs)
    for comb in combinations(vs, k):
        yield comb


def memory_expensive_random_choose(l, k, n=50000000):
    cvs = set()
    for elm in numpy.random.choice(list(l), (n, k), replace=True):
        elm = tuple(elm)
        if len({*elm}) != k:
            continue
        if elm in cvs:
            continue
        cvs.add(elm)
        yield elm
