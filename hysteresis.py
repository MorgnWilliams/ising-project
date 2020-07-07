import numpy as np
import matplotlib.pyplot as plt
from ising_lattice import IsingLattice


def M_as_function_of_H(H, T):
    # Computes magnetisation for each value of H in array
    N = 6
    M = np.zeros((len(H),))
    error = np.zeros((len(H),))

    for i in range(len(H)):
        if i == 0:
            # for the first H value, must do a number of steps to reach equilibrium
            steps_to_equilibrium = 30 * (N ** 2)
        else:
            # otherwise, value of H not much different from the previous value so only a few steps needed to equilibrate
            steps_to_equilibrium = 1 * (N ** 2)

        number_of_time_steps = 40 * (N ** 2)
        # calculate magnetisation and error for each value of H
        M[i], error[i] = lattice.compute_equilibrium_M(T, steps_to_equilibrium, number_of_time_steps)
        # Plot magnetisation per spin
        M[i] = np.divide(M[i], N ** 2)
        error[i] = np.divide(M[i], N ** 2)
    return M


def energy_as_function_of_H(H, T):
    # Computes energy for each value of H in array
    N = 6
    E = np.zeros((len(H),))
    error = np.zeros((len(H),))
    for i in range(len(H)):
        lattice._H = H[i]
        # generate a new lattice for each value of N
        if i == 0:
            # for the first H value, must do a number of steps to reach equilibrium
            steps_to_equilibrium = 30 * (N ** 2)
        else:
            # otherwise, value of H not much different from the previous value so only a few steps needed to equilibrate
            steps_to_equilibrium = 1 * (N ** 2)

        number_of_time_steps = 40 * (N ** 2)
        # calculate energy and error for each value of H
        E[i], error[i] = lattice.compute_equilibrium_energy(T, steps_to_equilibrium, number_of_time_steps)
        # Plot magnetisation per spin
        E[i] = np.divide(E[i], N ** 2)
        error[i] = np.divide(E[i], N ** 2)
    return E


# Plot the energy against H for 4 temperatures, three below and one at the critical temperature
T = np.array([0.5, 1.0, 1.5, 2.269])

for i in range(len(T)):
    # apply the magnetic field in the forward direction
    H_forward = np.linspace(-3, 3, 100)

    # initialise a new lattice for each value of temperature
    lattice = IsingLattice(6, 1, H_forward[0])

    E_forward = energy_as_function_of_H(H_forward, T[i])

    # reverse the direction of the application of the magnetic field
    H_backward = -H_forward
    E_backward = energy_as_function_of_H(H_backward, T[i])

    # append both forward and backward arrays for plotting
    E_total = np.append(E_forward,E_backward)
    H_total = np.append(H_forward, H_backward)
    plt.plot(H_total, E_total, '.', label = 'T = ' + str(T[i]) + r' $J/k_B$')

plt.title(r'Energy Against Applied Field for N = 6')
plt.xlabel(r'H / $1/\mu$')
plt.ylabel('Energy per Spin / J')
plt.legend()
plt.show()


# Plot the magnetisation against H for 4 temperatures, three below and one at the critical temperature
T = np.array([0.5, 1.0, 1.5, 2.269])

for i in range(len(T)):
    # apply the magnetic field in the forward direction
    H_forward = np.linspace(-3, 3, 100)

    # initialise a new lattice for each value of temperature
    lattice = IsingLattice(6, 1, H_forward[0])

    M_forward = M_as_function_of_H(H_forward, T[i])

    # reverse the direction of the application of the magnetic field
    H_backward = -H_forward
    M_backward = M_as_function_of_H(H_backward, T[i])

    # append both forward and backward arrays for plotting
    M_total = np.append(M_forward, M_backward)
    H_total = np.append(H_forward, H_backward)
    plt.plot(H_total, M_total, '.', label='T = ' + str(T[i]) + r' $J/k_B$')

plt.title(r'Magnetisation Against Applied Field for N = 6')
plt.xlabel(r'H / $1/\mu$')
plt.ylabel('Magnetisation per Spin')
plt.legend()
plt.show()
