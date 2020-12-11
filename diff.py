"""
Created on Aug 09 14:13:09 2019
@author: paul
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Importing data
data = np.loadtxt('20mT.txt', delimiter='\t')
size = len(data[:, 0])

# Converting to angular frequency
T_mod = data[:, 0]
omega = (2*np.pi) / (data[:, 0])
temp = (data[:, 1]) / 1000
omega = np.array(omega)
temp = np.array(temp)


# Fitting function (KVK theory)
def teta(x, d, G):
    # Bohr radius
    a_B = 1e-8
    # New designations
    a = - 2 / a_B
    k = - np.sqrt(x / (2 * d))
    b = k
    c = - np.sqrt(2 * x / d)

    # Fitting formula
    return 1 / ((16 * np.pi / a_B ** 3) * ((G / np.sqrt(d * x)) *
                (40 / a ** 6 + 6 / (a ** 4 * (a + b) ** 2) - 40 /
                ((a + c) ** 6) - 6 / ((a + b) ** 2 * (a + c) ** 4))) /
                np.sqrt(2))


# Fit (KVK theory)
popt, pcov = curve_fit(teta, omega, temp,
                       bounds=([1e-20, 1e15], [1e-15, 1e20]))

# Fit plot
omega_r = np.arange(0, max(omega), 0.1)
T_mod_r = np.arange(0, max(T_mod), 0.1)

# Figure plot
plt.plot(omega, temp, 'ro')
plt.plot(omega_r, teta(omega_r, *popt), 'b-')
plt.plot(omega_r, teta(omega_r, 9.99e-17, 3.5e16), 'g-')
plt.xlabel('omega (rad/s)')
plt.ylabel('teta (K)')
# plt.savefig('5mT.pdf')
plt.show()

# Fit parameters (D, G)
print('Fit parameters [D,G]:')
print(popt)

# Save fit to file
DataOut = np.column_stack((1 / (omega_r / (2 * np.pi)),
                           1000 * teta(omega_r, *popt)))
np.savetxt('diff_fit.dat', DataOut)
