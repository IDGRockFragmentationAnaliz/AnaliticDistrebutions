import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial, binom


def main():
    a = 0.5
    n_max = 30
    t = 1

    c = np.zeros((n_max+1, 1))

    for n in range(n_max+1):
        ms = np.arange(n+1)
        b1 = ((-1)**ms)*binom(n, ms)
        b2 = np.exp(-(a/2)*(ms - 1)*ms*np.log(2) - 2**(-a*(n - ms))*t)
        b = b1*b2
        c[n] = 2**n*np.sum(b)/(factorial(n)*(1-2**(-a))**n)


    plt.plot(c)
    plt.show()


if __name__ == '__main__':
    main()