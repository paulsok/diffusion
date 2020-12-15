"""
Created on Tue May 22 01:13:09 2018
@author: paul
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import csv

B = [None] * 200
rho = [None] * 200
BN_T = [None] * 14
time = np.arange(0, 1750, 125)


def func(x, S0, A, BN, b):
    Bz = 1.9  # mT
    ksi = 3
    BL = 0.2  # (mT) local field
    B1 = 4  # (mT) B_{1/2}
    x0 = - 1.3
    return (S0 * (Bz ** 2 / ((x + x0) ** 2 + Bz ** 2) + ((x + x0) ** 2 /
            ((x + x0) ** 2 + Bz ** 2)) * A / (1 + (np.sqrt((x + x0) ** 2 +
            Bz ** 2) - BN * Bz / np.sqrt((x + x0) ** 2 + Bz ** 2 +
            ksi * BL ** 2)) ** 2 / B1 ** 2)) + b)


for i in range(0, 14):
    data = np.loadtxt('/Users/paulsokolov/Downloads/Python/14.04.18/' +
                      repr(i + 1) + 'H.dat', skiprows=1, delimiter='\t')
    for j in range(0, 200):
        B[j] = data[j, 1]
        rho[j] = data[j, 6]

    B = np.array(B)
    rho = np.array(rho)
    popt, pcov = curve_fit(func, B, rho)
    BN_T[i] = popt[2]

    plt.plot(B, rho, 'ro')
    plt.plot(B, func(B, *popt), 'b-',
             label='S0=%5.3f, A=%5.3f, $B_N=%5.3f$ mT, b=%5.3f' % tuple(popt))
    plt.legend(loc=4)
    plt.xlabel('$B_x$ (mT)')
    plt.ylabel('rho')
    # plt.savefig(repr(i+1)+'H.png')
    plt.show()

with open('/Users/paulsokolov/Downloads/Python/sigma.dat', 'a') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(BN_T))
