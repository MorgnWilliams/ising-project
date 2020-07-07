import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def func(x, a, b, c):
    return a + b*(x**-(1/c))

N = np.array([6,7,8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20])
Tc_data = np.array([2.4256968916486565, 2.4012387116055955, 2.38649365398545, 2.3756893603020646, 2.3564256551160314, 2.347828208637039,  2.34704568197239, 2.3370919804181574, 2.330926924974846, 2.3266843121627407, 2.32706469192118, 2.314906655544282, 2.314335886])
Tc_errors = np.array([0.005022404371119432, 0.0024390206072329363,0.0065414812787397, 0.006999824190775233, 0.010070250299094667, 0.0032996552592925166, 0.0032196433100195, 0.008150357874463884,  0.0027742903099080424, 0.004329623855657852, 0.0120122484671728,  0.00460521371910873, 0.003372829])
popt, pcov = curve_fit(func, N, Tc_data, p0=[2.269, 1.0, 1.0])
perr = np.sqrt(np.diag(pcov))
plt.errorbar(N, Tc_data,yerr=Tc_errors, xerr=None, capsize= 5, ls='None', label = 'Data')
plt.title(r'Finite Size Scaling Plot to Determine $T_C(\infty)$')
plt.plot(N, func(N, *popt), label = r'$T_C(\infty) = 2.27 \pm 0.01$')
plt.xlabel('N')
plt.ylabel(r'$T_C(N)$ / K')
plt.legend()
plt.show()
print(popt)
print(perr)