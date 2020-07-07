import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from ising_lattice import IsingLattice


def func(x, a, b):
    # function to be fit using the curve_fit method
    return a * np.exp(-(1/b) * x)


def equilibrium_magnetisation_as_function_of_steps(N, T, max_tau):
    # calculate the magnetisation as a function of the number of steps after reaching equilibrium

    # steps to equilibrium as a function of temperature for each N is informed by part 1
    steps_to_equilibrium = 1000
    number_of_time_steps = max_tau

    # generate new lattice on which to compute autocovariance
    lattice = IsingLattice(N, 1, 0)

    # compute magnetisation per spin
    M = lattice.magnetisation_function_of_steps(T, steps_to_equilibrium + number_of_time_steps)
    # remove the steps prior to equilibration
    M = M[steps_to_equilibrium:]
    return M


def autocorrelation(x):
    # calculate auto correlation of a given function
    result = np.correlate(x, x, mode='full')
    result = result[len(x)-1:]
    return np.divide(result, result[0])


def auto_correlation_of_magnetisation(N, T, max_tau):
    # calculate the autocorrelation of the magnetisation and its error for a number of repeats
    number_of_repeats = 5
    auto_correlation_repeats = np.zeros((number_of_repeats, max_tau))

    for i in range(number_of_repeats):
        M = equilibrium_magnetisation_as_function_of_steps(N, T, max_tau)
        M = np.abs(M)
        # calculate autocovariance
        averageM = np.mean(M)
        auto_correlation_repeats[i] = autocorrelation(M - averageM)

    error = np.std(auto_correlation_repeats, axis=0)
    auto_correlation = np.mean(auto_correlation_repeats, axis=0)
    return auto_correlation, error


N = np.array([8, 12])
T = 1.5
tau = np.arange(0, 10002)
for i in range(len(N)):
    M = equilibrium_magnetisation_as_function_of_steps(N[i], T, len(tau))
    tau = tau[:10000]
    M = M[:10000]
    plt.plot(tau, M, '.')
    plt.title('Plot of Magnetisation as a Function of Time for N = ' + str(N[i]) + ', T = ' + str(T) + r' $J/k_B$')
    plt.xlabel('Time / Steps')
    plt.ylabel('Magnetisation Per Spin')
    plt.ylim([-1.1, 1.1])
    plt.show()

# Plot autocorrelation against steps for various N values on a semilog plot to observe divergence of tau_e at the
# critical point
T = 2.269
tau = np.arange(0, 10000)
colours = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red']
N = np.array([4, 8, 16])
tau_max = [5, 10, 15]

for i in range(len(N)):
    auto_correlation, error = auto_correlation_of_magnetisation(N[i], T, len(tau))
    auto_correlation = np.abs(auto_correlation)
    error = np.abs(error)
    # find tau_e by linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(tau[:tau_max[i]], np.log(auto_correlation[:tau_max[i]]))
    # plot the obtained data points with their errors
    plt.errorbar(tau[:tau_max[i]], auto_correlation[:tau_max[i]], yerr=error[:tau_max[i]], xerr=None,
                 capsize=5, label=str(N[i]), marker='.', ls='none', color=colours[i])
    # plot the obtained linear regression fit
    plt.plot(tau[:tau_max[i]], np.exp(intercept)*np.exp(tau[:tau_max[i]]*slope), color=colours[i])
plt.legend()
# set the scale of the y axis to log
plt.yscale('log', nonposy='clip')
plt.ylabel(r'a($\tau$)')
plt.xlabel(r'$\tau$ / Steps')
plt.title('Autocorrelation Against Lag Time for N = 4, 8, 16' + ', T = ' + str(T) + r' $J/k_B$')
plt.show()


# Plot autocorrelation against steps for each N value
tau = np.arange(0, 10000)
for i in range(len(N)):
    auto_correlation, error = auto_correlation_of_magnetisation(N[i], T, len(tau))
    plt.plot(tau,auto_correlation, '.', label='Data', markersize=3)
    plt.xlabel(r'$\tau$ / Steps')
    plt.ylabel(r'a($\tau$)')
    plt.title('Autocorrelation Against Lag Time for N = ' + str(N[i]) + ', T = ' + str(T) + r' $J/k_B$')

    # fit an expontential curve to determine tau_e and its error
    popt, pcov = curve_fit(func, tau, auto_correlation)
    perr = np.sqrt(np.diag(pcov))

    error = np.format_float_positional(perr[1], precision=2)
    tau_e = np.format_float_positional(popt[1], precision=2)

    # plot the curve fit over the data
    plt.plot(tau, func(tau, *popt), label=r'$\tau_e = $' + str(tau_e) + ' ' + r'$\pm$' + ' ' + str(error))
    plt.legend()
    plt.show()



