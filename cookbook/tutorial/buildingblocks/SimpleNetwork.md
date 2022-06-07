A Chemistry Solver from Scratch
================================

This tutorial walks you through the steps to create a new chemical kinetic rate
equation solver.  These solvers all utilize very similar structures for
calculating the right hand side and the Jacobian, and will export a
standardized API to those functions.  

While Dengo provides a number of pre-packaged rate equations, species, rate
coefficients and cooling functions, here we provide the necessary steps in
order to create new species, reaction rates, and a complete ODE solver from scratch.

```{tableofcontents}
```
