import itertools as it
import numpy as np

def gen_random_3(fineness=10):
    """

    For three variables, a Cartesian product enumeration with one
    condition is equivalent to geometric construction by an equilateral
    triangle to uniformly sample the composition space (create a uniform
    grid).

    """
    sqrt3 = 3**(1/2)
    det = 2*sqrt3
    xs = []
    linspace = [x / fineness for x in range(fineness + 1)]
    for d1 in linspace:
        for d2 in linspace:
            t1 = -sqrt3*(1-d1-2*d2) + sqrt3*(1-d1)
            t1 /= det
            t2 = sqrt3*(1-d1-2*d2) + sqrt3*(1-d1)
            t2 /= det

            assert round(t1, 3) == round(d2, 3)
            assert round(t2, 3) == round(1 - d1 - d2, 3)

            if 0 <= round(t1,3) <= 1 and 0 <= round(t2,3) <= 1:
                ycoord = (1 - d1 - d2)*sqrt3
                t3 = (ycoord * 2 / sqrt3) / 2
                assert round(t3, 3) == round(1 - d1 - d2, 3)
                xs.append((d1, d2, round(1 - d1 - d2, 3)))
    return xs

def gen_random_n(n,fineness=10):
    combs = it.product(range(fineness + 1),repeat=n)
    keep = []
    for comb in combs:
        if sum(comb) == fineness:
            keep.append(tuple(x/fineness for x in comb))
    return keep

def gen_random_n_fast(n, K=None):
    """

    Generates random samples n-dimensional data from the flat Dirichlet distribution.

    Input: 
      n: the dimension of the composition space desired
      K: the number of random samples desired
    Output:
      samples: K-by-n array containing K samples
    """ 
    if K is None:
        K = 10**n
    # overdraw, since you will have to discard many samples
    x1 = np.random.random(10*n*K)
    x2 = np.random.random(10*n*K)
    bl = np.exp(-x1) > x2
    samples = x1[bl][:n*K]
    samples = samples.reshape((K,n)) 
    norm = samples.sum(axis=1)
    return samples / norm[:, np.newaxis]
