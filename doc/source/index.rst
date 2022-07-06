.. Dengo documentation master file, created by
   sphinx-quickstart on Thu May  9 10:57:47 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Dengo: An Engineer for Chemistry Solvers
========================================

Hi there!  Welcome to Dengo.  Dengo is a Python system for symbolically
describing a system of chemical species, chemical kinetic rate equations,
cooling functions, and then producing from these a set of numerical kernels to
calculate the "right hand side" and Jacobian of this system.  These two
numerical kernels can then be linked into one of several ODE solvers.

Dengo is best thought of as a way to create a chemistry and cooling solver that
you can stick into a codebase.  Rather than trying to tie together the four
separate parts -- species, reactions, rate coefficients and the solver itself
-- into a single whole, Dengo allows you to construct each item individually
and join them at the final step before inserting them into a simulation code.

Contents:

.. toctree::
   :maxdepth: 2

   tutorial

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
