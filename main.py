import numpy as np
import matplotlib.pyplot as plt



def main():
    # n = 10
    # ks =  np.arange(n+1)
    #
    # B = np.zeros((len(ks), 1))
    # for i in range(len(ks)):
    #     B[i] = np.exp(get_lnB(i, n, 0.5))
    #
    # plt.plot(ks, B)
    # plt.show()

    n_max = 40
    a = 0.25
    ns = np.arange(n_max+1)

    t = 1000

    N = np.zeros((len(ns), 1))
    N2 = np.zeros((len(ns), 1))

    for n in range(n_max+1):
        N[n] = get_Nn(n, a, t=t)
        N2[n] = get_Nn2(n, a, t=t)

    #

    #for k in range(n+1):
    #    lnB[k] = get_lnB(k, n, a)

    lam = np.log(a*t)/(np.log(2)*a)
    print(lam)
    plt.plot(ns, N)
    plt.plot(ns, N2)
    plt.plot([lam, lam], [0, np.max(N)])
    plt.show()


def get_Nn2(n, a, t=0):
    lnB = get_lnB(n, n, a)
    N = np.exp(lnB - 2**(-a*n)*t)
    return N

def get_Nn(n, a, t=0):
    N = 0
    for k in range(n+1):
        lnB = get_lnB(k, n, a)
        B = (-1)**(n-k)*np.exp(lnB - 2**(-a*k)*t + n)
        N = N + B
    return N

def get_lnB(k, n, a):
    lnB1 = lnQP2(a, k)
    lnB2 = lnQP2(a, n - k)
    lnE = -(a / 2) * (n - k - 1) * (n - k) * np.log(2)
    return lnE + lnB1 + lnB2

def lnQP2(a, k, max_j=100):
    # Инициализация суммы
    result = -a*k**2/2*np.log(2)
    return result

def lnQP(a, k, max_j=100):
    # Инициализация суммы
    result = 0.0
    # Вычисление суммы
    for j in range(1, max_j + 1):
        term = (1 - 2 ** (-a * k * j)) / (1 - 2 ** (-a * j)) * (2 ** (-a * j) / j)
        result += term
    return result


if __name__ == '__main__':
    main()