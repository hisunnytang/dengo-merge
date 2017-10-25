Tutorial: A Chemistry Solver from Scratch
=========================================

This tutorial walks you through the steps to create a new chemical kinetic rate
equation solver.  These solvers all utilize very similar structures for
calculating the right hand side and the Jacobian, and will export a
standardized API to those functions.  

While Dengo provides a number of pre-packaged rate equations, species, rate
coefficients and cooling functions, here we provide the necessary steps in
order to create new species and rates from scratch.

Defining a Species
------------------

We start by defining individual species.  This can be done inside a python
module of your choosing, which we will run at the end to output our new
network.  Species are defined by a small number of attributes:

 * Name (which will be used in cooling functions and internally to the solver)
 * Number: Mostly unused except when handling ions.
 * Atomic weight (in integer AMU)
 * Number of free electrons that is contributes

This information is used when calculating things like the contribution of a
species to the number density.

.. warning:: At the present time, Dengo does not support variable gammas for
   species.  This is a planned future incorporation.

To create a new species, you can both create the species object *and* register
it in the global dictionary like so:

.. code-block:: python

   from dengo.reaction_classes import Species

   HI = Species('HI', 1.0, 1.0, 0.0)
   HII = Species("HII", 1.0, 1.0, 1.0)
   de = Species("de", 1.0, 1.0, 0.0)

We now have three symbolic "species" objects for hydrogen, ionized hydrogen,
and electrons.  Note that Dengo will happily create ions for species defined in
the CHIANTI database.

Creating Reactions
------------------



Specifying Reaction Rate Coefficients
-------------------------------------

Cooling Functions
-----------------

Creating a Network
------------------
