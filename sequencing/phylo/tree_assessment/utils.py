import random
import numpy
import scipy
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


def ordered_memory_expensive_random_choose(l, k, n=50000000):
    cvs = set()
    for elm in numpy.random.choice(list(l),(n,k),replace=True):
        elm = tuple(elm)
        if len({*elm}) != k:
            continue
        if elm in cvs:
            continue
        cvs.add(elm)
        yield elm


def memory_expensive_random_choose(l, k, n=50000000):
    l = list(l)
    n = min(int(scipy.special.comb(len(l), k, repetition=True)), n)
    cvs = set()
    for elm in numpy.random.choice(l,(n,k),replace=True):
        elm = tuple(elm)
        if len({*elm}) != k:
            continue
        elm = frozenset(elm)
        if elm in cvs:
            continue
        cvs.add(elm)
        yield elm