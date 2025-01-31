import numpy as np
import matplotlib.pyplot as plt
from pandas import Float32Dtype
from scipy.special import factorial


from sage.all import *
from sage.combinat.q_analogues import gaussian_binomial
from sage.combinat.q_analogues import q_factorial
from sage.combinat.q_analogues import q_pochhammer


def main():
    a = 0.2
    n_max = 40
    t = 500
    q = 2 ** (-a)
    c = np.zeros((n_max+1, 1))

    for n in range(n_max+1):
        ks = np.arange(n+1)
        binom = np.array([gaussian_binomial(Integer(n), Integer(k), q=RealNumber(q)) for k in range(n+1)])
        b1 = ((-1)**(n - ks))*binom
        b2 = np.exp((1/2)*((n - ks) - 1)*(n - ks)*np.log(q) - (q**ks)*t)
        b = b1*b2
        qf = q_factorial(n, q=RealNumber(q))*((1-q)**n)
        c[n] = np.sum(b)/qf

    s = np.sum(c)
    print(s)
    plt.plot(c)
    plt.show()


if __name__ == '__main__':
    main()