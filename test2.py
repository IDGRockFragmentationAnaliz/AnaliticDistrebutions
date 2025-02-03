import numpy as np
import matplotlib.pyplot as plt
from pandas import Float32Dtype
from scipy.special import factorial


from sage.all import *
from sage.combinat.q_analogues import gaussian_binomial
from sage.combinat.q_analogues import q_factorial
from sage.combinat.q_analogues import q_pochhammer



def main():
    a = 1
    n_max = 20
    t = 100
    q = 2 ** (-a)
    c = np.zeros((n_max+1, 1))

    for n in range(n_max+1):
        c[n] = get_cn(n, t, q)

    s = np.sum(c)
    print(s)

    #mu = np.log(a*t)/(a*np.log(2))

    mu = get_d(1, 1, q)
    print(mu)
    exit()
    plt.plot([mu, mu],[0, np.max(c)],color="r")
    plt.plot(c,color="black")
    plt.show()


def get_mu(q, tau):
    sumK1 = 0
    sumK2 = 0
    for k in range(100):
        q_poch = q_pochhammer(Integer(k), RealNumber(q), RealNumber(q))
        dsum = np.exp(-(q**k)*tau)/q_poch
        sumK1 += dsum
        sumK2 += k*dsum
    print(sumK1, sumK2)


def get_u(i, q, m_max=100):
    result = 0
    for m in range(m_max):
        qp = q_pochhammer(Integer(m), RealNumber(q), RealNumber(q))
        result += (-1)**m * q**((m-1)*m/2)/qp * (m**i)
    return result

def get_d(k, i, q):
    bin = binomial_coefficients(i)
    result = 0
    for s in range(i):
        result += (bin[(i, s)]*(k**(i - s))*get_u(s, q))
    return result


def get_cn(n, t, q):
    ks = np.arange(n + 1)
    bin1 = np.array([gaussian_binomial(Integer(n), Integer(k), q=RealNumber(q)) for k in range(n + 1)])
    b1 = ((-1) ** (n - ks)) * bin1
    b2 = np.exp((1 / 2) * ((n - ks) - 1) * (n - ks) * np.log(q) - (q ** ks) * t)
    b = b1 * b2
    qf = q_factorial(n, q=RealNumber(q)) * ((1 - q) ** n)
    return np.sum(b) / qf


if __name__ == '__main__':
    main()