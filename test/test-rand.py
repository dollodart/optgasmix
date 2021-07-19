from optgasmix.rand import gen_random_3,\
        gen_random_n, gen_random_n_fast
from time import time

def test_gen_random_3(display = False):
    t0 = time()
    xs = gen_random_3()
    dt = time() - t0

    if display:
        x1, x2, x3 = zip(*xs)
        x1 = sorted(x1)
        x2 = sorted(x2)
        x3 = sorted(x3)

        for x, y, z in zip(x1, x2, x3):
            print(x, y, z)
    return dt

def test_gen_random_n(display = False):

    t0 = time()
    xs = gen_random_n(4,10)
    dt = time() - t0

    if display:
        x1, x2, x3, x4 = zip(*xs)
        x1 = sorted(x1)
        x2 = sorted(x2)
        x3 = sorted(x3)
        x4 = sorted(x4)
        for t in zip(x1,x2,x3,x4):
            print(' '.join(str(x) for x in t))
    return dt


def test_gen_random_n_fast(display=False):

    t0 = time()
    xs = gen_random_n_fast(4,10**4)
    dt = time() - t0

    if display: 
        import matplotlib.pyplot as plt
        # inspect histograms for uniform random values
        plt.hist(xs[:,0])
        plt.hist(xs[:,1])
        plt.hist(xs[:,2])
        plt.show() 
        # visually inspect tetrahedron (expensive)
        #from mpl_toolkits.mplot3d import Axes3D
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #ax.scatter(xs[:,0], xs[:,1], xs[:,2])
        #plt.show()
    return dt

if __name__ == '__main__':
    dt = test_gen_random_3()
    print(dt)
    dt1 = test_gen_random_n()
    dt2 = test_gen_random_n_fast()
    print(dt1, dt2)
