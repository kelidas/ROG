import math
import itertools
import numpy as np

def rog_lengths(n, nvar):
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

    Returns
    -------
    lengths : array dtype float
        pairwise distances among points
    counts : array of dtype int
        number of distances of the same type

    Examples
    --------
    >>> n = 2
    >>> nvar = 3
    >>> rog_lengths(n, nvar)
    (array([ 0.5 , 0.70710678 , 0.8660254 ]), \
     array([12, 12,  4]))
    '''
    # prepare array of possible rank number differences
    d = np.arange(n)
    lengths = np.array([], dtype=int)
    counts = np.array([], dtype=int)
    # prepare iterator for difference vectors delta_t
    deltas = itertools.combinations_with_replacement(d, nvar)
    next(deltas) # skip the first null vector (0,...,0)
    # loop over difference vectors delta_t
    for delta_t in deltas:
        delta_t = np.array(delta_t) # convert list to array
        # sum of squared differences
        length_t = np.sum(delta_t**2)

        # number n_t^n    eq. ((*@\ref{eq:FFD-coeff_n_i^n}@*))
        ntn = np.prod((n - delta_t))
        # number n_t^p    eq. ((*@\ref{eq:FFD-coeff_n_i^p}@*))
        ntp = math.factorial(nvar)
        for k in set(delta_t):
            # frequencies of differences in the vector delta_t
            ntp /= math.factorial(np.sum((delta_t - k)==0))
        # number n_t^d eq. ((*@\ref{eq:FFD-coeff_n_i^d}@*))
        ntd = 2 ** (np.sum(delta_t > 0) - 1)
        # number n_t eq. ((*@\ref{eq:nt}@*))
        nt = ntn * ntp * ntd

        lengths = np.append(lengths, length_t)
        counts = np.append(counts, nt)
    # calculate real lengths in unit hypercube eq. ((*@\ref{eq:L_ij:delta_ij}@*))
    lengths = np.sqrt(lengths) / float(n)
    return lengths, counts

if __name__ == '__main__':
    n = 2
    nvar = 30
    print(rog_lengths(n, nvar))
