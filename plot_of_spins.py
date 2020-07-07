import numpy as np
import matplotlib.pyplot as plt
from ising_lattice import IsingLattice

f = plt.figure(figsize=(15, 15), dpi=80);
T = np.array([1.5, 2.269, 3.0])
for i in range(len(T)):
    lattice = IsingLattice(32, 1, 0)
    for j in range(10000):
        lattice.perform_metropolis(T[i])
    lattice.plot_configuration(f, T[i], 32, i+1)
plt.show()
