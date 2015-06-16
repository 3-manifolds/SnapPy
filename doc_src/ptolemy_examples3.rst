Step-by-step examples: Part 3
=============================

.. _ptolemy-example-using-magma-sage:

Computing solutions with magma or sage vs retrieving solutions
--------------------------------------------------------------

So far, we querried the database for solutions to a Ptolemy variety::

    >>> p = Manifold("m011").ptolemy_variety(2)
    >>> p.retrieve_solutions()
    Trying to retrieve solutions from http://ptolemy.unhyperbolic.org/data/pgl2/OrientableCuspedCensus/03_tetrahedra/m011__sl2_c0.magma_out ...
    Parsing...
    [PtolemyCoordinates(
         {'c_0011_0': 1,
          'c_0011_1': -1,
          'c_0011_2': -1,
          'c_0101_0': -1,
          'c_0101_1': Mod(x^2 + x, x^3 + 2*x^2 + x + 1),
	  ...
          's_3_2': 1},
         is_numerical = False, ...)]

We can use ``compute_solutions`` instead of ``retrieve_solutions`` to actually compute the solutions ourselves (for example, for a non-census triangulation not in the database). Currently, we support two engines: 

* `sage <http://www.sagemath.org/>`_ (which is free, but can only solve a fairly limited number of Ptolemy varieties)
* `magma <http://magma.maths.usyd.edu.au/magma/>`_

If you are inside sage::

    >>> p = Manifold("m011").ptolemy_variety(2)
    >>> p.compute_solutions(engine = 'sage')
    [PtolemyCoordinates(
         {'c_0011_0': 1,
          'c_0011_1': -1,
          'c_0011_2': -1,
          'c_0101_0': -1,
          'c_0101_1': Mod(x^2 + x, x^3 + 2*x^2 + x + 1),
	  ...
          's_3_2': 1},
         is_numerical = False, ...)]

If you have magma installed::

    >>> p = Manifold("m011").ptolemy_variety(2)
    >>> p.compute_solutions(engine = 'magma', verbose = True)
    Writing to file: /tmp/tmppNSc8S/m011__sl2_c0.magma
    Magma's output in: /tmp/tmppNSc8S/m011__sl2_c0.magma_out
    Command: ulimit -m 732421; echo | magma "/tmp/tmppNSc8S/m011__sl2_c0.magma" > "/tmp/tmppNSc8S/m011__sl2_c0.magma_out"
    Starting magma...
    magma finished.
    Parsing magma result...
    [PtolemyCoordinates(
         {'c_0011_0': 1,
          'c_0011_1': -1,
          'c_0011_2': -1,
          'c_0101_0': -1,
          'c_0101_1': Mod(x^2 + x, x^3 + 2*x^2 + x + 1),
	  ...
          's_3_2': 1},
         is_numerical = False, ...)]

To get an idea of what Ptolemy varieties magma can still handle, have a look at the `database <http://ptolemy.unhyperbolic.org/html/summary.html>`_: for *N* = 2, the computations up to 12 tetrahedra only took 

**Remark:** The magma engine is not expected to work under windows. It will also fail if magma is not installed or the magma executable cannot be found. The ptolemy module creates a temporary file (``m011__sl2_c0`` here) and also gives the command it tried to run to process the file through magma. If you believe that magma is installed correctly on your system but encounter an error, you can try to run the command yourself to understand better what is going on. Feel free to report a bug (to enischte at gmail dot com) including the temporary files (``m011__sl2_c0`` and ``m011__sl2_c0.out`` here) and any other error messages.

**Remark:** If no engine is specified, it is assumed to be sage when used inside sage and magma instead.

.. _ptolemy-example-complex-volume:

Computing the complex volume
----------------------------

Similar to ``volume_numerical``, we can compute the complex volume (volume + i Chern-Simons) for all representations (that are :ref:`generically decorated <ptolemy-generically-decorated>`).

Here is an example computing the solutions to the Ptolemy variety ourselves::

    >>> Manifold("m011").ptolemy_variety(2,'all').compute_solutions().complex_volume_numerical()
    [[[-4.30211422042248 E-16 + 0.725471193740844*I,
       -0.942707362776931 + 0.459731436553693*I,
       0.942707362776931 + 0.459731436553693*I]],
     [[4.64255370258293 E-15 + 0.680993020093457*I,
       3.94215909915729 E-15 + 0.312682687518267*I,
       -2.78183391239608 - 0.496837853805869*I,
       2.78183391239608 - 0.496837853805869*I]]]

And here the same example retrieving solutions from the database::

    >>> Manifold("m011").ptolemy_variety(2,'all').retrieve_solutions().complex_volume_numerical()
    Trying to retrieve solutions from http://ptolemy.unhyperbolic.org/data/pgl2/OrientableCuspedCensus/03_tetrahedra/m011__sl2_c0.magma_out ...
    Parsing...
    Trying to retrieve solutions from http://ptolemy.unhyperbolic.org/data/pgl2/OrientableCuspedCensus/03_tetrahedra/m011__sl2_c1.magma_out ...
    Parsing...
    [[[-4.30211422042248 E-16 + 0.725471193740844*I,
    ...
       2.78183391239608 - 0.496837853805869*I]]]

.. _ptolemy-detailed-example-matrices:

Computing the matrices for a different presentation
---------------------------------------------------

The ``fundamental_group`` method of a SnapPy triangulation can yield different presentations by supplying  optional arguments such as ``simplify_presentation`` and ``minimize_number_of_generators``. If we have a word in one of these presentations and want to evaluate its image under the representation, we need to supply the presentation as follows::

    >>> M=Manifold("m003")
    >>> sol = M.ptolemy_variety(2).retrieve_solutions()[0]
    >>> G = M.fundamental_group(simplify_presentation = False)
    >>> sol.evaluate_word('a', G)

Again, we can check that the representation actually assigns the identity to all relators:

    >>> G
    Generators:
       a,b,c
    Relators:
       BCaC
       AbCbA
    >>> sol.evaluate_word('AbCbA', G)
    [[Mod(1, x^2 - x - 1), 0], [0, Mod(1, x^2 - x - 1)]]
    >>> for relator in G.relators():
    ...     print sol.evaluate_word(relator, G)
    [[Mod(1, x^2 - x - 1), 0], [0, Mod(1, x^2 - x - 1)]]
    [[Mod(1, x^2 - x - 1), 0], [0, Mod(1, x^2 - x - 1)]]
  

.. _ptolemy-example-boundary-holonomy:

Computing the images of the peripheral curves for a representation
------------------------------------------------------------------

The object returned by ``fundamental_group`` also contains words for the peripheral curves of a manifold. We can compute the corresponding matrices::

    >>> M = Manifold("m003")
    >>> G = M.fundamental_group()
    >>> sol = M.ptolemy_variety(2,1).retrieve_solutions()[0]
    >>> for i, cusp_curves in enumerate(G.peripheral_curves()):
    ...     print "Cusp %d:" % i
    ...     for cusp_curve in cusp_curves:
    ...         print sol.evaluate_word(cusp_curve, G)
    Cusp 0:
    [[Mod(2*x - 3, x^2 + x + 1), Mod(2*x, x^2 + x + 1)], [Mod(6, x^2 + x + 1), Mod(-2*x + 1, x^2 + x + 1)]]
    [[Mod(-2*x - 5, x^2 + x + 1), Mod(-2, x^2 + x + 1)], [Mod(6*x + 6, x^2 + x + 1), Mod(2*x + 3, x^2 + x + 1)]]    

**Remark:** For each cusp, we can conjugate these matrices into *P* since the representation is :ref:`boundary-unipotent <ptolemy-boundary-unipotent>`. We might implement a method returning a matrix in *P* for the longitude and meridian of a cusp in the future (simply by finding loops corresponding to a longitude and meridian as path of short edges in the truncated simplex in Figure 17 of [GGZ2012]_).

.. _ptolemy-example-find-witness:

Finding a witness point for a positively dimensional compoent of the Ptolemy variety
------------------------------------------------------------------------------------

We already saw an :ref:`example of a positively dimensional component <ptolemy-non-zero-dim-comp>`. By flattening, we obtain a list of all the components of the Ptolemy varieties for all obstruction classes::

    >>> M=Manifold("m371")
    >>> sols = M.ptolemy_variety(2,'all').retrieve_solutions().flatten()

We can now just look at the positively dimensional ones::

    >>> one_dim_sols = [ sol for sol in sols if sol.dimension > 0]
    >>> len(one_dim_sols)
    1
    >>> one_dim_sols
    [[ PtolemyCoordinates(
           {'c_0011_0': 1,
            'c_0011_1': -1,
            'c_0011_2': -1,
            'c_0011_3': Mod(-x - 1, x^2 + x + 2),
	    ...,
            's_3_4': 1},
           is_numerical = False, ...) 
       (witnesses for NonZeroDimensionalComponent(dimension = 1, free_variables = ['c_0110_2'])) ]]

We see that we have one such component and that each component is actually itself a list of witness points.

**Remark:** Witness points are a fairly new feature and not all files in the database have been updated yet to contain them. You might instead just see ``[NonZeroDimensionalComponent(dimension = 1)]``.

**Remark:** The ptolemy module also reports the `free variables` for the positively dimensional components. We can set these variables to random values and generically will obtain a new witness point.

We can access the witness point(s) for each component just by iteration::

    >>> for component in one_dim_sols:
    ...     print "Component:"
    ...     for witness in component:
    ...         print "    Witness:"
    ...         print "        Volumes:", witness.volume_numerical()
    Component:
        Witness:
            Volumes: [0.E-38, 0.E-38]

The different volumes in a line correspond to different Galois conjugates of the same `algebraic` witness point. 

.. _ptolemy-example-non-zero-dim-rep:

Finding non-zero dimensional families of boundary-unipotent representations
---------------------------------------------------------------------------

We now revisit the :ref:`1-dimensional component of the Ptolemy variety<ptolemy-non-zero-dim-comp>` and answer the question whether this yields a 1-dimensional family of representations or not. We pick a :ref:`witness point <ptolemy-example-find-witness>` for the component and check the :ref:`matrices for the peripheral curves <ptolemy-example-boundary-holonomy>`::

    >>> M = Manifold("m371")
    >>> G = M.fundamental_group()
    >>> sols = M.ptolemy_variety(2,'all').retrieve_solutions().flatten()
    >>> components = [ sol for sol in sols if sol.dimension > 0]
    >>> for component in components:
    ...     print "Component of dimension %d" % component.dimension
    ...     for witness in component:
    ...         for i, cusp_curves in enumerate(G.peripheral_curves()):
    ...             print "    Cusp %d:" % i
    ...             for cusp_curve in cusp_curves:
    ...                 print "        ", witness.evaluate_word(cusp_curve, G)
    Component of dimension 1
        Cusp 0:
            [[Mod(1, x^2 + x + 2), 0], [0, Mod(1, x^2 + x + 2)]]
            [[Mod(1, x^2 + x + 2), 0], [0, Mod(1, x^2 + x + 2)]]

We see that the matrices are trivial, thus this 1-dimensional component corresponds to a 1-dimensional family of :ref:`decorations <ptolemy-generically-decorated>` of the same (up to Galois conjugacy) representation. The corresponding family of representation is 0-dimensional.

Let us try another manifold, ``m410``:

    >>> M = Manifold("m410")
    >>> G = M.fundamental_group()
    >>> sols = M.ptolemy_variety(2,'all').retrieve_solutions().flatten()
    >>> components = [ sol for sol in sols if sol.dimension > 0]
    >>> for component in components:
    ...     print "Component of dimension %d" % component.dimension
    ...     for witness in component:
    ...         for i, cusp_curves in enumerate(G.peripheral_curves()):
    ...             print "    Cusp %d:" % i
    ...             for cusp_curve in cusp_curves:
    ...                 print "       ", witness.evaluate_word(cusp_curve, G)
    Component of dimension 1
        Cusp 0:
            [[Mod(1, x^2 + 2), 0], [0, Mod(1, x^2 + 2)]]
            [[Mod(1, x^2 + 2), Mod(x, x^2 + 2)], [0, Mod(1, x^2 + 2)]]
    Component of dimension 1
        Cusp 0:
            [[Mod(1, x^2 + 7), 0], [0, Mod(1, x^2 + 7)]]
            [[Mod(1, x^2 + 7), 0], [0, Mod(1, x^2 + 7)]]

It has two 1-dimensional components, and for the first one, we see that the matrices are non-trivial, so this corresponds indeed to a 1-dimensional family of representations.

**Remark:** The witness points are chosen so that they are not at the intersection of two positively dimensional components. This is for the following reason: it could happen that there is a 1-dimensional family of representations which contains points where the boundary holonomy becomes trivial. This yields a representation where the above matrices are trivial yet it is part of a 1-dimensional family of boundary-unipotent representations. In the ptolemy variety, however, this means that two non-zero dimensional components (one corresponding to a family of decorations, the other to a family of representations) intersect.

Representations that are the same as PSL(2, **C**)-representations
------------------------------------------------------------------

Let us compare the volumes of ``m009`` and ``m159``:

    >>> Manifold("m009").ptolemy_variety(2,'all').retrieve_solutions().volume_numerical()
    [[],
     [],
     [],
     [[2.66674478344907, -2.66674478344907, 2.66674478344907, -2.66674478344907]]]
    >>> Manifold("m159").ptolemy_variety(2,'all').retrieve_solutions().volume_numerical()
    [[[0.E-38, 0.E-37, 0.E-37],
     [-2.02988321281931, 2.02988321281931, -2.02988321281931, 2.02988321281931]],
    [[0.698544082784440, -0.698544082784440, 3.82168758617998, -3.82168758617998],
     [0.E-37, 0.E-37, 0.E-37, 0.E-37]]]

In both cases, some volumes appear twice (2.66674... for ``m009`` and 2.02988... for ``m159``). In the case of ``m009``, these two volumes correspond to the same PSL(2, **C**)-representation and in case ``m159`` to two different boundary-unipotent SL(2, **C**)-representations that are the same as PSL(2, **C**)-representations (see :ref:`ptolemy-psl-multiplicity`). We can get the  multiplicity by calling ``degree_to_shapes``::

    >>> Manifold("m009").ptolemy_variety(2).degree_to_shapes()
    2
    >>> Manifold("m159").ptolemy_variety(2).degree_to_shapes()
    1

When we convert the Ptolemy coordinates to shapes/cross ratios for ``m009``, we also see that we see the same shape assignment appears twice (at least numerically)::

    >>> Manifold("m009").ptolemy_variety(2,'all').retrieve_solutions().numerical().cross_ratios()	
    [[],
     [],
     [],
     [[CrossRatios(
           {'z_0000_0': 0.500000000000000 + 1.32287565553230*I,
            'z_0000_1': 0.375000000000000 + 0.330718913883074*I,
            'z_0000_2': 0.500000000000000 + 1.32287565553230*I,
            'zp_0000_0': 0.250000000000000 + 0.661437827766148*I,
    	    ...},
           is_numerical = True, ...),
       CrossRatios(
           {'z_0000_0': 0.500000000000000 - 1.32287565553230*I,
            'z_0000_1': 0.375000000000000 - 0.330718913883074*I,
            'z_0000_2': 0.500000000000000 - 1.32287565553230*I,
    	    ...}
           is_numerical = True, ...),
       CrossRatios(
           {'z_0000_0': 0.500000000000000 + 1.32287565553230*I,
            'z_0000_1': 0.375000000000000 + 0.330718913883074*I,
            'z_0000_2': 0.500000000000000 + 1.32287565553230*I,
            ...},
           is_numerical = True, ...),
       CrossRatios(
           {'z_0000_0': 0.500000000000000 - 1.32287565553230*I,
            'z_0000_1': 0.375000000000000 - 0.330718913883074*I,
            'z_0000_2': 0.500000000000000 - 1.32287565553230*I,
  	    ...},
           is_numerical = True, ...)]]]
    


**Remark:** The tables at `ptolemy.unhyperbolic.org <http://ptolemy.unhyperbolic.org/html/summary.html>`_ use the cross ratios to list representations that are the same as PSL(2, **C**)-representation only once.

