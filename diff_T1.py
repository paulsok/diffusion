"""
Created on Aug 09 14:13:09 2019
@author: paul
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Importing data
data = np.loadtxt('5mT.txt', delimiter='\t')

# Converting to angular frequency
T_mod = data[:, 0]  # (s)
omega = 2 * np.pi / T_mod  # (rad/s)
temp = (data[:, 1]) / 1000  # (K)


# Fitting function (KVK theory)
def teta(x, D, G):
    # Bohr radius
    a_B = 10*1e-9  # (nm)
    # T1 nuclear relaxation time
    T1 = 10  # (s)
    # New designations
    a = -2 / a_B
    eta = ((1 / (2 * np.sqrt(D * T1))) * ((np.sqrt(1 +
           (x * T1) ** 2) - 1) / (1 + (x * T1) ** 2) ** 0.25))
    ksi = ((1 / (2 * np.sqrt(D * T1))) * ((np.sqrt(1 + (x * T1) ** 2) + 1) /
                                          (1 + (x * T1) ** 2) ** 0.25))
    b = - (ksi + 0 * eta)
    # Fitting formula
    return - 1 / (((12 * np.pi * G) / (a_B ** 3 * D)) *
                  ((14 * a ** 4 + 20 * a ** 3 * b + 15 * a ** 2 * b ** 2 + 6 *
                    a * b ** 3 + b ** 4) / (a ** 5 * (a + b) ** 6)))
    # return 1 / ((3 * G * (ksi ** 2 - eta ** 2)) /
    #             (32 * a_B * D * (ksi ** 2 + eta ** 2) ** 2))


def teta_r(x, D, G):
    # Bohr radius
    a_B = 10 * 1e-9  # (nm)
    # T1 nuclear relaxation time
    T1 = 10  # (s)
    # New designations
    eta = (1 / (2 * np.sqrt(D * T1))) * ((np.sqrt(1 + (x * T1) ** 2) - 1) /
                                         (1 + (x * T1) ** 2) ** 0.25)
    ksi = (1 / (2 * np.sqrt(D * T1))) * ((np.sqrt(1 + (x * T1) ** 2) + 1) /
                                         (1 + (x * T1) ** 2) ** 0.25)
    # Fitting formula
    return 1 / ((3 * G * (ksi ** 2 - eta ** 2)) / (32 * a_B * D *
                (ksi ** 2 + eta ** 2) ** 2))


# Fit (KVK theory)
popt, pcov = curve_fit(teta, omega, temp, bounds=([1e-25, 1], [2e-14, 1e21]))

omega_r = np.arange(min(omega), max(omega), 0.01)

# Figure plot
plt.plot(omega, temp, 'ro')
plt.plot(omega_r, teta(omega_r, *popt), 'b-')
# plt.plot(omega_r, teta_r(omega_r, 1e-14,1e1), 'k-')
plt.xlabel('omega (rad/s)')
plt.ylabel('theta (K)')
plt.grid(True)
plt.show()
# plt.savefig('5mT.png')

# Fit parameters (D,G)
print('Fit parameters [D,G]:')
print(popt)

# Save fit to file
# DataOut_automatic = np.column_stack(((2 * np.pi) /
#                     omega_r, 1000 * teta(omega_r, *popt)))
# np.savetxt('diff_fit_automatic.dat', DataOut_automatic )
