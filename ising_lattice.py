import numpy as np
import matplotlib.pyplot as plt


class IsingLattice:

    def __init__(self, N, J, H):
        # constructor method for the Ising Lattice
        self._N = N
        self._J = J
        self._H = H
        # initialise the lattice spins to a random array on instantiation of the lattice object
        self.initial_spins()

    def initial_spins(self):
        # initialises the lattice spins to a random array
        # create N by N random 2D array of 0, 1
        self._spins = np.random.randint(2, size=(self._N, self._N))
        # replace all instances of 0 with -1, faster than looping
        self._spins = np.where(self._spins == 0, -1, self._spins)

    def perform_metropolis(self, T):
        # perform one iteration of the Metropolis algorithm
        for i in range(self._N):
            for j in range(self._N):
                # loop through N**2 points
                # generate random lattice site (a, b) for each loop iteration
                a = np.random.randint(0, self._N)
                b = np.random.randint(0, self._N)
                # compute the change in energy on flipping the spin
                dE = self.compute_dE(a, b)
                # generate a random number such that the spin may be flipped with the Boltzmann probability
                p = np.random.uniform(0,1)
                if dE < 0:
                    # energetically favourable, flip
                    self._spins[a, b] *= -1
                elif np.exp(-dE/T) > p:
                    # energetically unfavourable but flip with probability exp(-dE/T)
                    self._spins[a, b] *= -1
                else:
                    # energetically unfavourable reject the flip with probability 1-exp(-dE/T)
                    pass

    def compute_dE(self, i, j):
        # compute the energy of a given spin flip
        N = self._N
        H = self._H
        # apply the periodic boundary conditions with modulo operators to bring all 'neighbours' outside the lattice
        # inside
        E = -self._J*(self._spins[(i+1) % N, j] + self._spins[i, (j+1) % N] +
                      self._spins[(i-1) % N, j] + self._spins[i, (j-1) % N])
        dE = -2*E*self._spins[i, j] + 2*H*self._spins[i, j]
        return dE

    def check_uniform(self):
        # check whether the spins are in one of the two ground states
        s = self._spins
        # flatten the array to a 1D C-style array
        s = s.ravel()
        if np.count_nonzero(s == s[0]) == len(s):
            return True
        else:
            return False

    def compute_equilibrium_M(self, T, steps_to_equilibrium, number_of_time_steps):
        # function calculates the equilibrium magnetisation for a given temperature using 100 bins of statistically
        # independent configurations
        self.initial_spins()  # remove if computing hysteresis

        # perform a number of steps to reach equilibrium
        for i in range(0, steps_to_equilibrium):
            self.perform_metropolis(T)

        # system now in equilibrium, can now compute magnetisation
        # use 100 bins for measurements
        number_of_bins = 100
        steps_per_bin = number_of_time_steps//number_of_bins
        # initialise magnetisation array
        M = np.zeros((number_of_bins, steps_per_bin), dtype=int)

        for i in range(number_of_bins):
            for j in range(steps_per_bin):
                M[i][j] = np.sum(self._spins)
                self.perform_metropolis(T)
        # find the mean of each bin value
        M_bin_averages = np.mean(M, axis=1)
        # error on the final mean value of M
        error = np.std(M_bin_averages)
        # average over the bin values
        M = np.mean(M_bin_averages)
        return M, error

    def magnetisation_function_of_steps(self, T, max_steps):
        # returns magnetisation per spin as a function of number of steps
        self.initial_spins()
        M = np.zeros((max_steps,), dtype=int)
        for i in range(0, max_steps-1):
            self.perform_metropolis(T)
            M[i] = np.sum(self._spins)
        M = np.divide(M, self._N**2)
        return M

    def calculate_energy(self):
        # calculate the energy of a given orientation of spins
        s = self._spins
        N = self._N
        H = self._H
        energy = 0
        for i in range(len(s)):
            for j in range(len(s)):
                S = s[i, j]
                nb = s[(i + 1) % N, j] + s[i, (j + 1) % N] + s[(i - 1) % N, j] + s[i, (j - 1) % N]
                energy += -nb * S
        return (energy / 4) - H * np.sum(s)

    def compute_equilibrium_energy(self, T, steps_to_equilibrium, number_of_time_steps):
        # calculate the energy of an equilibrium configuration at a given T, H
        self.initial_spins()    # remove if computing hysteresis

        # perform a number of steps to reach equilibrium
        for i in range(0, steps_to_equilibrium):
            self.perform_metropolis(T)

        # system now in equilibrium, can now compute magnetisation
        # use 100 bins for measurements
        number_of_bins = 100
        steps_per_bin = number_of_time_steps//number_of_bins
        energy = np.zeros((number_of_bins, steps_per_bin), dtype=int)
        for i in range(number_of_bins):
            for j in range(steps_per_bin):
                energy[i][j] = self.calculate_energy()
                self.perform_metropolis(T)

        # average within the bins
        energy_bin_averages = np.mean(energy, axis=1)
        error = np.std(energy_bin_averages)
        # average over the bin values
        energy = np.mean(energy_bin_averages)
        return energy, error

    def compute_heat_capacity(self, T, steps_to_equilibrium, number_of_time_steps):
        # compute the heat capacity in equilibrium for a given T, H using the fluctuation-dissipation theorem result
        energy_bin_averages, error = self.compute_equilibrium_energy(T, steps_to_equilibrium, number_of_time_steps)
        C = (error**2)/((T**2))
        return C

    def plot_configuration(self, f, T, N, n):
        # Plot the spins for a given T, N
        spins = self._spins
        X, Y = np.meshgrid(range(N), range(N))
        sp = f.add_subplot(1, 3, n)
        plt.setp(sp.get_yticklabels(), visible=False)
        plt.setp(sp.get_xticklabels(), visible=False)
        plt.pcolormesh(X, Y, spins, cmap=plt.cm.RdBu)
        plt.title('Temperature = ' + str(T) + ' K')
        plt.axis('tight')
    plt.show()




