import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from ising_lattice import IsingLattice


def func(x, a, b, c):
    # Lorentzian to fit to heat capacity data
    return a * (b / ((x - c) ** 2 + b ** 2))


def heat_capacity_as_function_of_T(N, T, steps_to_equilibrium, number_of_time_steps):
    # computes heat capacity for a given temperature
    lattice = IsingLattice(N, 1, 0)
    C = np.mean(lattice.compute_heat_capacity(T, steps_to_equilibrium, number_of_time_steps))
    return C


N = [4, 5, 6]
heat_capacity_as_function_of_T = np.vectorize(heat_capacity_as_function_of_T)
T = np.linspace(1.5, 3.0, 50)

# lists in which to store values of the transition temperature and its errors
Tc = []
Tc_errors = []
number_of_repeats = 1

colours = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:purple']
for i in range(len(N)):
    C = []
    Tc_repeats = []
    for j in range(number_of_repeats):
        # calculate the heat capacity as function of T, and append this to array
        C_iteration = heat_capacity_as_function_of_T(N[i], T, 1000, 2000)
        C.append(C_iteration)
        # fit a Lorentzian curve to this heat capacity
        popt, pcov = curve_fit(func, T, C_iteration, p0=[1.0, 1.0, 2.269])
        # fitting of Lorentzian enables determination of Tc
        Tc_repeats.append(popt[2])
    # find the average and error of the repeats
    Tc.append(np.mean(Tc_repeats))
    Tc_errors.append(np.std(Tc_repeats))

    C = np.mean(C, axis=0)
    plt.plot(T, C, '.', label='N = ' + str(N[i]), color=colours[i])
    popt, pcov = curve_fit(func, T, C, p0 =[1.0, 1.0, 2.269])
    plt.plot(T, func(T, *popt), color=colours[i])

plt.title('Heat Capacity as a Function of T')
plt.xlabel(r'T / $\frac{J}{k_B}$')
plt.ylabel(r'$\frac{C}{k_B}$')
plt.legend()
plt.show()
