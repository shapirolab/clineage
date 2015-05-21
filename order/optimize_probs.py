from numpy.linalg import matrix_power
from numpy import matrix
from scipy.stats import binom
from collections import defaultdict
import numpy as np

# multiple up/ multiple dw implementation
def markovian_matrix(fus, fds, n=100):
    """
    :param fus: list of fu probability functions from +1 to +n
    :param fds: list of fd probability functions from -1 to -n
    :param n:
    :return:
    """
    assert len(fus) <= n/2
    assert len(fds) <= n/2

    def fs(n):
        return 1-sum([fu(n) for fu in fus])-sum([fd(n) for fd in fds])

    l = []
    for i in xrange(len(fds)):  # first rows with fds squeezing
        l.append(
            [fd(i) for fd in reversed(fds[:i])] +
            [sum([fd(i) for fd in fds[i:]])+fs(i)] +
            [fu(i) for fu in fus] + [0]*(n-i-len(fus)-1)
        )
    for i in xrange(len(fds), n-len(fus)):  # intermediate rows without squeezing
        l.append([0]*(i-len(fds)) +
                 [fd(i) for fd in reversed(fds)] +
                 [fs(i)] +
                 [fu(i)for fu in fus] +
                 [0]*(n-i-len(fus)-1))

    for fui, i in enumerate(xrange(n-len(fus)+1, n+1)):  # last rows with fus squeezing
        l.append([0]*(n-len(fds)-(len(fus)-fui)) + [fd(i) for fd in reversed(fds)] +
                 [fs(n) + sum([fu(i) for fu in fus[len(fus)-fui-1:]])] +
                 [fu(i) for fu in fus[:len(fus)-fui-1]])
    return matrix(l)


def get_probs_n(mat, n):
    return [mat[n, i] for i in xrange(mat.shape[0])]


def dyn_mat_model(fus, fds, ms_len, cycles, n=100):
    return get_probs_n(matrix_power(markovian_matrix(fus, fds, n=n), cycles), ms_len)