def gen_random_3(fineness=10):
    """

    For three variables, a Cartesian product enumeration with one condition is equivalent to
    geometric construction by an equilateral triangle to uniformly sample the composition space (create a uniform grid). 

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

if __name__ == '__main__':
    xs = gen_random_3()
    x1, x2, x3 = zip(*xs)
    x1 = sorted(x1)
    x2 = sorted(x2)
    x3 = sorted(x3)

    for x, y, z in zip(x1, x2, x3):
        print(x, y, z)

#xs = []
#for d1 in linspace:
#    for d2 in linspace:
#        if 1 - d1 - d2 >=0 :
#            for d3 in linspace:
#                if 1 - d1 - d2 - d3 >= 0:
#                    xs.append((d1, d2, d3, round( 1 - d1 - d2 - d3, 3)))
#
#x1, x2, x3, x4 = zip(*xs)
#x1 = sorted(x1)
#x2 = sorted(x2)
#x3 = sorted(x3)
#x4 = sorted(x4)
#for t in zip(x1,x2,x3,x4):
#    print(' '.join(str(x) for x in t))
