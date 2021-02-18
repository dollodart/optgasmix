import numpy as np
def mixrule(mu, mm, x):
    """
    Equation (1.4-15) and (9.3-17) of BSL.
    """
    mumix = 0
    n = len(x)
    for a in range(n):
        phib = []
        for b in range(n):
            mmr = mm[a] / mm[b]
            mur = mu[a] / mu[b]
            res = (1/8)**(1/2)*(1 + mmr)**(-1/2)*(1 + mur**(1/2)*mmr**(1/4))**2
            phib.append(res)
        mumix += x[a]*mu[a]/sum(phib[i] * x[i] for i in range(n))
    return mumix

def weighted_geometric_mean(x, w):
    """
    Weighted geometric mean. 
    """
    x = np.array(x)
    w = np.array(w)
    return np.exp(sum(w*np.log(x))/sum(w))

def weighted_arithmetic_mean(x, w):
    """
    Weighted arithmetic mean.
    """
    w = np.array(w)
    x = np.array(x)
    return sum(w*x)

def weighted_harmonic_mean(x, w):
    """
    Weighted harmonic mean.
    """
    x = np.array(x)
    w = np.array(w)
    return 1./(sum(w/x)/sum(w))
