import numpy as np
import matplotlib.pyplot as plt
from ising_lattice import IsingLattice

T = np.linspace(0.01, 5, 30)

# compute magnetisation
N = np.array([4, 8])
for i in range(len(N)):
    # generate a new lattice for each value of N
    lattice = IsingLattice(N[i], 1, 0)
    steps_to_equilibrium = 10 * (N[i]**2)
    number_of_time_steps = 100 * (N[i]**2)
    M = np.zeros((len(T),))
    error = np.zeros((len(T),))
    for j in range(len(T)):
        # calculate magnetisation and error for each temperature
        M[j], error[j] = lattice.compute_equilibrium_M(T[j], steps_to_equilibrium, number_of_time_steps)
    M = np.abs(M)
    # Plot magnetisation per spin
    M = np.divide(M, N[i]**2)
    error = np.divide(M, N[i]**2)

    plt.errorbar(T, M, yerr=error, xerr=None, ls='', label='N = ' + str(N[i]), fmt='.', capsize=2)

plt.legend()
plt.ylabel('Magnetisation per Spin')
plt.xlabel(r'T / $\frac{J}{k_B}$')
plt.title('Magnetisation Against Temperature for N = 4, 8')
plt.show()
