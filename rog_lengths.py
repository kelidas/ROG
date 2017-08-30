import math
import itertools
import numpy as np
from functools import reduce
import time


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print ('%r  %2.2f ms' % (method.__name__, (te - ts) * 1000))
        return result

    return timed


def rog_lengths_numpy(n, nvar, periodic=False):
    '''
    Evaluate pairwise distances among points forming
    regular orthogonal grid in a hypercube.

    Parameters
    ----------
    n : int
        number of equidistant points
        along an individual dimension
    nvar : int
        number of input random variables
        (dimension of a hypercube)
    periodic : bool, optional
        evaluate lengths in a periodically extended space.
        Default is False.

    Returns
    -------
    lengths : array dtype float
        pairwise distances among points
    counts : array of dtype int
        number of distances of the same type

    Examples
    --------
    >>> n = 3
    >>> nvar = 2
    >>> rog_lengths_numpy(n, nvar, periodic=False)
    (array([ 0.33333333,  0.66666667,  0.47140452,  0.74535599,  0.94280904]), \
array([12,  6,  8,  8,  2]))
     >>> rog_lengths_numpy(n, nvar, periodic=True)
     (array([ 0.33333333,  0.33333333,  0.47140452,  0.47140452,  0.47140452]), \
array([12,  6,  8,  8,  2]))
    '''
    lengths = []
    counts = []
    # prepare iterator for difference vectors delta_t
    deltas = itertools.combinations_with_replacement(range(n), nvar)
    next(deltas)  # skip the first null vector (0,...,0)
    # loop over difference vectors delta_t
    for delta_t in deltas:
        delta_t = np.array(delta_t)  # convert list to array
        # sum of squared differences
        if periodic:
            # update delta_t for periodic space ((*@\ref{eq:delta_bar}@*))
            h = (delta_t - n // 2) > 0  # Heaviside function
            delta_t_pae = np.abs(h * n - delta_t)
            length_t = np.sum(delta_t_pae**2)
        else:
            length_t = np.sum(delta_t**2)

        # number n_t^n    eq. ((*@\ref{eq:FFD-coeff_n_i^n}@*))
        ntn = np.prod((n - delta_t))
        # number n_t^p    eq. ((*@\ref{eq:FFD-coeff_n_i^p}@*))
        ntp = math.factorial(nvar)
        for k in set(delta_t):
            # frequencies of differences in the vector delta_t
            ntp //= math.factorial(np.count_nonzero((delta_t - k) == 0))
        # number n_t^d eq. ((*@\ref{eq:FFD-coeff_n_i^d}@*))
        ntd = 2 ** (np.count_nonzero(delta_t > 0) - 1)
        # number n_t eq. ((*@\ref{eq:nt}@*))
        nt = ntn * ntp * ntd

        lengths.append(length_t)
        counts.append(nt)
    # calculate real lengths in unit hypercube eq. ((*@\ref{eq:L_ij:delta_ij}@*)) or ((*@\ref{eq:L^bar_ij:delta_ij}@*))
    lengths = np.sqrt(lengths) / float(n)
    return lengths, np.asarray(counts)


def rog_lengths(n, nvar, periodic=False):
    '''
    Evaluate pairwise distances among points forming
    regular orthogonal grid in a hypercube.

    Parameters
    ----------
    n : int
        number of equidistant points
        along an individual dimension
    nvar : int
        number of input random variables
        (dimension of a hypercube)
    periodic : bool, optional
        evaluate lengths in a periodically extended space.
        Default is False.

    Returns
    -------
    lengths : list of floats
        pairwise distances among points
    counts : list of ints
        number of distances of the same type

    Examples
    --------
    >>> n = 3
    >>> nvar = 2
    >>> rog_lengths(n, nvar, periodic=False)
    ([0.3333333333333333, 0.6666666666666666, 0.47140452079103173, \
0.7453559924999299, 0.9428090415820635], [12, 6, 8, 8, 2])
     >>> rog_lengths(n, nvar, periodic=True)
     ([0.3333333333333333, 0.3333333333333333, 0.47140452079103173, \
0.47140452079103173, 0.47140452079103173], [12, 6, 8, 8, 2])
    '''
    lengths = []
    counts = []
    # prepare iterator for difference vectors delta_t
    deltas = itertools.combinations_with_replacement(range(n), nvar)
    next(deltas)  # skip the first null vector (0,...,0)
    # loop over difference vectors delta_t
    for delta_t in deltas:
        # sum of squared differences
        if periodic:
            # update delta_t for periodic space ((*@\ref{eq:delta_bar}@*))
            h = [(dt - n // 2) > 0 for dt in delta_t]
            delta_t_pae = [abs(hi * n - dt) for hi, dt in zip(h, delta_t)]
            length_t = sum([dt ** 2 for dt in delta_t_pae])
        else:
            length_t = sum([dt ** 2 for dt in delta_t])

        # number n_t^n    eq. ((*@\ref{eq:FFD-coeff_n_i^n}@*))
        ntn = reduce(lambda x, y: x * y, [n - dt for dt in delta_t])
        # number n_t^p    eq. ((*@\ref{eq:FFD-coeff_n_i^p}@*))
        ntp = math.factorial(nvar)
        for k in set(delta_t):
            # frequencies of differences in the vector delta_t
            ntp //= math.factorial(sum([(dt - k) == 0 for dt in delta_t]))
        # number n_t^d eq. ((*@\ref{eq:FFD-coeff_n_i^d}@*))
        ntd = 2 ** (sum([dt > 0 for dt in delta_t]) - 1)
        # number n_t eq. ((*@\ref{eq:nt}@*))
        nt = ntn * ntp * ntd

        lengths.append(length_t)
        counts.append(nt)
    # calculate real lengths in unit hypercube eq. ((*@\ref{eq:L_ij:delta_ij}@*)) or ((*@\ref{eq:L^bar_ij:delta_ij}@*))
    lengths = [l ** 0.5 / float(n) for l in lengths]
    return lengths, counts


if __name__ == '__main__':
    n = 20
    nvar = 3

    # TEST equality of results for both alternatives
    l, c = timeit(rog_lengths_numpy)(n, nvar, periodic=False)
    l2, c2 = timeit(rog_lengths)(n, nvar, periodic=False)
    print(np.allclose(l, l2))
    print(np.allclose(c, c2))

    l, c = timeit(rog_lengths_numpy)(n, nvar, periodic=True)
    l2, c2 = timeit(rog_lengths)(n, nvar, periodic=True)
    print(np.allclose(l, l2))
    print(np.allclose(c, c2))

    # documentation tests
    import doctest
    doctest.testmod()
