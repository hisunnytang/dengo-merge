# Compare Dengo-Enzo Integration with Enzo-Grackle

[$\texttt{Grackle}$](https://github.com/grackle-project/grackle) is a chemistry and cooling library for astrophysical simulations and modelling. It allows for non-equilibrium primordial chemistry calculations.
It contains innate support to a wide range of simulation softwares.
In this section, we aim to compare our $\texttt{Enzo}$ integration with $\texttt{Grackle}$. The simulation, as we will discuss further in this chapter, is initialized with cosmological initial conditions, and is evolved with the 9-species primordial chemistry model. The simulation is centered on a $10^6 M_\odot$ dark matter halo at around redshift 16.5, the simulation is halted when the central density reaches $\sim 10^{-10} \mathrm{g ~cm^{-3}}$.

A brief introduction to the processes involved is discussed at first, and is followed by an outline of how to integrate $\texttt{Dengo}$ with $\texttt{Enzo}$ Finally, with the help of [$\texttt{yt}$](https://yt-project.org/doc), we will compare the sliceplots, radial profiles, cooling/heating rates of these simulations.

```{tableofcontents}
```
