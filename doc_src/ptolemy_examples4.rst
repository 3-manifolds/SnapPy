Step-by-step examples: Part 4
=============================

**Advanced topics. Still under construction.**

.. _ptolemy-example-structure-of-solution:

The structure of an exact solution
----------------------------------

A solution to the Ptolemy variety (here ``sol``) is an object of type ``PtolemyCoordinates`` which is a subclass of a python dictionary. It assigns a value to each Ptolemy coordinate c\ :sub:`...` (the s\ :sub:`...` are related to the obstruction class which is trivial here). We can get the value assigned to a particular Ptolemy coordinate as follows::

    >>> sol = Manifold("m003").ptolemy_variety(2).retrieve_solutions()[0]
    >>> sol['c_0101_1']
    Mod(-x, x^2 - x - 1)

This is a pari object of type `POLMOD`. It means that the solution is in the number field with defining polynomial *x*\ :sup:`2`\ - *x* - 1 and that it is equal to - *x* where *x* is a root in the defining polynomial.

We can get to the parts of a pari `POLMOD` object by::

    >>> sol['c_0101_1'].lift()
    -x
    >>> sol['c_0101_1'].mod()
    x^2 - x - 1

These objects also support the field operations, here we are computing a cross ratio::

    >>> (sol['c_0101_0'] * sol['c_1010_0']) / (sol['c_1001_0'] * sol['c_0110_0'])
    Mod(x, x^2 - x - 1)

**Remark:** We would prefer to represent these using types that are better suited when using sage such as NumberField and NumberFieldElement. In the future, we might jettison support of the ptolemy module outside of sage and use sage's native types.

Rational Univariate Representation
----------------------------------

We have some limited support for the `Rational Univariate Representation`, for example::

    >>> sols = Manifold("m069").ptolemy_variety(3).retrieve_solutions()
    >>> sols[4]
    PtolemyCoordinates(
        { 'c_0201_0': ( Mod(-426088934700737884313191344*x^24 + 4110489425474123899213651272*x^23 ...
                                 - 67064980598091504185767190,
                            x^25 ... - 196124*x^4 + 14010*x^3 - 1560*x^2 + 72*x - 1)
                  ) / ( Mod(875895332415105303646551573*x^24 - 8450034810535061601312104296*x^23 ...
                                 + 137871639973525691285247446, 
                            x^25 ... - 196124*x^4 + 14010*x^3 - 1560*x^2 + 72*x - 1) ),
          ...
        },
        is_numerical = False, ...)

This means that ``c_0201_0`` is assigned the fraction of those two ``POLDMOD`` objects.

**Remark:** This is still under development.

TODO
----

* Give an overview of the classes and the methods.
* Explain how the magma file are generated. Say that magma computes the RadicalDecomposition. Ptolemy module can process the lexicographic Groebner basis of a prime ideal.


