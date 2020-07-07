import numpy as np
import matplotlib.pyplot as plt
from ising_lattice import IsingLattice

# initialise a figure in which to plot magnetisation as a function of number of steps for different N values
f = plt.figure(figsize=(15, 15), dpi=80, tight_layout=True)

T = 1.5
number_of_steps = np.arange(1, 100)
N = np.array([4, 8, 12])
for i in range(len(N)):
    # generate a new lattice object for each value of N
    lattice = IsingLattice(N[i], 1.0, 0)
    # calculate the magnetisation as a function of the number of steps
    M = lattice.magnetisation_function_of_steps(T, len(number_of_steps))
    f.add_subplot(3, 3, i+1)
    plt.plot(number_of_steps, M, '.')
    plt.title('N = ' + str(N[i]))
    plt.xlabel('Steps')
    plt.ylabel('Magnetisation per Spin')
plt.savefig('magnetisation.png')
plt.show()