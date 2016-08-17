from numpy import zeros, matrix


def build_markovian_matrix(steps, n=50, squeeze=False):
    mat = zeros((n, n))
    for step, func in steps.items():
        for i in range(n):
            val = func(i)
            if 0 <= i + step < n:
                mat[i, i + step] = val
            elif squeeze:
                ind = min(n-1, max(0, i+step))
                mat[i, ind] += val
    return mat
