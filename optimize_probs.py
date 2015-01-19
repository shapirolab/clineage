from numpy.linalg import matrix_power
from numpy import matrix

def markovian_matrix(fu, fd, n=100):
    def fs(n):
        return 1-fu(n)-fd(n)
    return matrix([[fd(0)+fs(0),fu(0)] + [0]*(n-2)] + 
                  [[0]*(i-1) + [fd(i), fs(i), fu(i)] + [0]*(n-i-2)
                    for i in xrange(1,n-1)] + 
                  [[0]*(n-2) + [fd(n), fs(n)+fu(n)]])

def get_probs_n(mat, n):
    return [mat[n,i] for i in xrange(mat.shape[0])]


def dyn_mat_model(fu, fd, ms_len, cycles):
    return get_probs_n(matrix_power(markovian_matrix(fu, fd), cycles),ms_len)
