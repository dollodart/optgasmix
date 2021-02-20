import itertools as it
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

import numpy as np
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

if __name__ == '__main__':
    from time import time
    t0 = time()
    xs = gen_random_3()
    print(time() - t0)
    x1, x2, x3 = zip(*xs)
    x1 = sorted(x1)
    x2 = sorted(x2)
    x3 = sorted(x3)

    for x, y, z in zip(x1, x2, x3):
        pass
        #print(x, y, z)

    t0 = time()
    xs = gen_random_n(4,10)
    print(time() - t0)
    x1, x2, x3, x4 = zip(*xs)
    x1 = sorted(x1)
    x2 = sorted(x2)
    x3 = sorted(x3)
    x4 = sorted(x4)
    for t in zip(x1,x2,x3,x4):
        pass
        #print(' '.join(str(x) for x in t))


    t0 = time()
    xs = gen_random_n_fast(4,10**4)
    print(xs[0])
    print(time() - t0)
    #import matplotlib.pyplot as plt
    # visually inspect tetrahedron (expensive)
    #from mpl_toolkits.mplot3d import Axes3D
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(xs[:,0], xs[:,1], xs[:,2])
    #plt.show()
    # inspect histograms for uniform random values
    #plt.hist(xs[:,0])
    #plt.hist(xs[:,1])
    #plt.hist(xs[:,2])
    #plt.show()
