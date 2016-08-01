from numpy.linalg import matrix_power
from numpy import matrix


def build_markovian_matrix(fus, fs, fds, n):
    assert len(fus) <= n / 2
    assert len(fds) <= n / 2

    l = []
    for i in range(len(fds)):  # first rows with fds squeezing
        l.append(
            [fd(i) for fd in reversed(fds[:i])] +
            [sum([fd(i) for fd in fds[i:]]) + fs(i)] +
            [fu(i) for fu in fus] + [0] * (n - i - len(fus) - 1)
        )
    for i in range(len(fds), n - len(fus)):  # intermediate rows without squeezing
        l.append([0] * (i - len(fds)) +
                 [fd(i) for fd in reversed(fds)] +
                 [fs(i)] +
                 [fu(i) for fu in fus] +
                 [0] * (n - i - len(fus) - 1))

    for fui, i in enumerate(range(n - len(fus) + 1, n + 1)):  # last rows with fus squeezing
        l.append([0] * (n - len(fds) - (len(fus) - fui)) + [fd(i) for fd in reversed(fds)] +
                 [fs(n) + sum([fu(i) for fu in fus[len(fus) - fui - 1:]])] +
                 [fu(i) for fu in fus[:len(fus) - fui - 1]])
    return matrix(l)


# multiple up/ multiple dw implementation
def stochastic_markovian_matrix(fus, fds, n=100):
    """
    :param fus: list of fu probability functions from +1 to +n
    :param fds: list of fd probability functions from -1 to -n
    :param n:
    :return:
    """
    def fs(n):
        return 1-sum([fu(n) for fu in fus])-sum([fd(n) for fd in fds])

    return build_markovian_matrix(fus, fs, fds, n)


# multiple up/ multiple dw implementation
def non_stochastic_markovian_matrix(fus, fds, n=100):
    fs = fds[0]
    fds = fds[1:]
    return build_markovian_matrix(fus, fs, fds, n)


def get_probs_n(mat, n):
    return [mat[n, i] for i in range(mat.shape[0])]


def dyn_mat_model(fus, fds, ms_len, cycles, n=100):
    return get_probs_n(matrix_power(stochastic_markovian_matrix(fus, fds, n=n), cycles), ms_len)

def dyn_mat_fs_model(fus, fds, ms_len, cycles, n=100):
    return get_probs_n(matrix_power(non_stochastic_markovian_matrix(fus, fds, n=n), cycles), ms_len)