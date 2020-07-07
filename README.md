# The Ising Model of Ferromagnetism
An investigation of the Ising Model of Ferromagnetism in 2D using the Metropolis Algorithm is presented.

* ising_lattice.py  contains the 2D Ising Lattice class
* equilibration.py  investigates the Monte-Carlo time for equilibration to a given temperature of a lattice of dimensions N x N
* autocorrelation.py  investigates the statistical independence of configurations as a function of the lag time between configurations, this is used to compute appropriate erros on physical quantities
* magnetisation.py  equipped with the ability to measure physical quantities and their errors, we investigate the variation of the magnetisation with temperature
* heat_capacity.py  using the fluctuation-dissipation theorem we determine the heat capacity as a function of temperature
* finite_size_scaling.py  investigates the finite-size scaling of the system
* hysteresis.py applying a magnetic field, and cycling its values at different temperatures, the hysteretic behaviour of the system below the transition temperature is observed

