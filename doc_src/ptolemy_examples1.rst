Step-by-step examples: Part 1
=============================

.. _ptolemy-example-basic:

The Ptolemy variety for SL(*N*, **C**)
--------------------------------------

Given a SnapPy triangulation, we obtain the reduced Ptolemy variety to find
SL(2, **C**)-representations as follows::

  >>> M=Manifold("m003")
  >>> M.ptolemy_variety(N = 2)
  Ptolemy Variety for m003, N = 2
    c_0011_0 * c_0101_0 + c_0011_0^2 - c_0101_0^2
    c_0011_0 * c_0101_0 + c_0011_0^2 - c_0101_0^2
    - 1 + c_0011_0

The result of ``M.ptolemy_variety(2)`` is an object of type ``PtolemyVariety``.

**Remark:** The exact formatting of the output might change between SnapPy versions and sage.

**Remark:** The first two equations are the two Ptolemy relations 
for the two tetrahedra in ``m003``. The last equation :ref:`reduces <ptolemy-reduced-variety>` the Ptolemy variety.

Similar, we can obtain the Ptolemy variety for higher *N*, say SL(3, **C**)::

  >>> M=Manifold("m004")
  >>> M.ptolemy_variety(3)
  Ptolemy Variety for m003, N = 3
        c_0012_0 * c_1101_0 + c_0012_1 * c_0111_0 - c_0102_0 * c_1011_0
        c_0012_1 * c_1110_0 - c_0102_0 * c_0111_0 + c_0102_1 * c_1011_0
        c_0012_0 * c_0111_0 + c_0102_0 * c_1101_0 - c_0102_1 * c_1110_0
        c_0012_0 * c_1110_0 + c_0012_1 * c_1011_0 - c_0102_1 * c_1101_0
        - c_0012_0 * c_0111_0 - c_0012_1 * c_1101_0 + c_0102_1 * c_1011_0
        - c_0012_0 * c_1110_0 - c_0102_0 * c_1011_0 + c_0102_1 * c_0111_0
        - c_0012_1 * c_0111_0 + c_0102_0 * c_1110_0 - c_0102_1 * c_1101_0
        - c_0012_0 * c_1011_0 - c_0012_1 * c_1110_0 + c_0102_0 * c_1101_0
        - 1 + c_0012_0
        - 1 + c_0111_0

**Remark:** Similarly, we obtain four Ptolemy relations for each of the two tetrahedra in ``m004`` corresponding to the four subsimplices of a tetrahedron we get for *N*\ =3 (see Figure 2 of [GTZ2011]_).

Using auto-completion
---------------------

Let us assign a Ptolemy variety to ``p`` and then type ``p.``::

    >>> p=Manifold("m003").ptolemy_variety(2)
    >>> p.

If we are in SnapPy, sage or ipython, we can now hit the tab-key and see a list of attributes and methods available for a Ptolemy variety::

    >>> p.
    p.canonical_representative           p.filename_base                      p.to_magma
    p.compute_decomposition              p.path_to_file                       p.to_magma_file
    p.compute_solutions                  p.py_eval_section                    p.variables
    p.degree_to_shapes                   p.py_eval_variable_dict              p.variables_with_non_zero_condition
    p.equations                          p.retrieve_decomposition             
    p.equations_with_non_zero_condition  p.retrieve_solutions                 

We can get further help by using the ``?``::

    >>> p.filename_base?
    ...
    Definition: p.filename_base(self)
    Docstring:
    Preferred filename base for writing out this Ptolemy variety
    ...

This is a general mechanism and works for all objects in SnapPy, sage or ipython.

.. _ptolemy-example-retrieve-exact-solutions:

Retrieving exact solutions from the database
--------------------------------------------

Given a Ptolemy variety, we can access the database at `ptolemy.unhyperbolic.org <http://ptolemy.unhyperbolic.org/>`_ to retrieve solutions for it with ``retrieve_solutions`` (if this is not working, please check your Internet connection)::

    >>> p=Manifold("m003").ptolemy_variety(2)
    >>> sols=p.retrieve_solutions()
    Trying to retrieve solutions from http://ptolemy.unhyperbolic.org/data/pgl2/OrientableCuspedCensus/02_tetrahedra/m003__sl2_c0.magma_out ...
    Parsing..
    >>> sols
    [PtolemyCoordinates(
         {'c_0011_0': 1,
          'c_0011_1': -1,
          'c_0101_0': Mod(x, x^2 - x - 1),
          'c_0101_1': Mod(-x, x^2 - x - 1),
          'c_0110_0': Mod(-x, x^2 - x - 1),
          'c_0110_1': Mod(x, x^2 - x - 1),
          'c_1001_0': -1,
          'c_1001_1': 1,
          'c_1010_0': Mod(x, x^2 - x - 1),
          'c_1010_1': Mod(-x, x^2 - x - 1),
          'c_1100_0': 1,
          'c_1100_1': -1,
          's_0_0': 1,
          's_0_1': 1,
          's_1_0': 1,
          's_1_1': 1,
          's_2_0': 1,
          's_2_1': 1,
          's_3_0': 1,
          's_3_1': 1},
         is_numerical = False, ...)]
    
The result is a list of solutions (up to Galois conjugation), here the list contains only one solution. Let us pick that one::

    >>> len(sols)
    1
    >>> sol = sols[0]
    PtolemyCoordinates(
        {'c_0011_0': 1,
         'c_0011_1': -1,
         'c_0101_0': Mod(x, x^2 - x - 1),
         'c_0101_1': Mod(-x, x^2 - x - 1),
         'c_0110_0': Mod(-x, x^2 - x - 1),
         'c_0110_1': Mod(x, x^2 - x - 1),
         'c_1001_0': -1,
         'c_1001_1': 1,
         'c_1010_0': Mod(x, x^2 - x - 1),
         'c_1010_1': Mod(-x, x^2 - x - 1),
         'c_1100_0': 1,
         'c_1100_1': -1,
         's_0_0': 1,
         's_0_1': 1,
         's_1_0': 1,
         's_1_1': 1,
         's_2_0': 1,
         's_2_1': 1,
         's_3_0': 1,
         's_3_1': 1},
        is_numerical = False, ...)

As we can see, a solution assigns a value to each Ptolemy coordinate c\ :sub:`...`\ . It is of type ``PtolemyCoordinates`` (a subclass of python's ``dict``) and more details are discussed in :ref:`a later example <ptolemy-example-structure-of-solution>`.

**Remark:** We can give the additional argument ``verbose=False`` to suppress the messages about the database access::

    >>> sols = Manifold("m003").ptolemy_variety(2).retrieve_solutions(verbose=False)

.. _ptolemy-example-matrices:

Compute the matrices for a representation
-----------------------------------------

**Remark:** Requires SnapPy 2.3 or later.

Given a solution as above, we can take a word in the fundamental group and get its image under the representation using ``evaluate_word``. Here, we do it for the two generators::

    >>> M = Manifold("m003")
    >>> sol = M.ptolemy_variety(2).retrieve_solutions()[0]
    >>> sol.evaluate_word('a')
    [[0, Mod(1, x^2 - x - 1)], [Mod(-1, x^2 - x - 1), Mod(-x, x^2 - x - 1)]]
    >>> sol.evaluate_word('b')
    [[Mod(x, x^2 - x - 1), Mod(x, x^2 - x - 1)],
    [Mod(-x, x^2 - x - 1), Mod(-1, x^2 - x - 1)]]

By default, this word is with respect to the presentation of the fundamental group that SnapPy computes when given no further arguments. Thus, we expect the identity matrix when we evaluate a relator (for PSL(*N*, **C**) the diagonal element will be an *N*-th root of unity)::

    >>> M.fundamental_group()
    Generators:
       a,b
    Relators:
       abAAbabbb
    >>> sol.evaluate_word('abAAbabbb')
    [[Mod(1, x^2 - x - 1), 0], [0, Mod(1, x^2 - x - 1)]]
    
We revisit computing the matrices :ref:`here <ptolemy-detailed-example-matrices>` to explain how to use a different presentation of the fundamental group.

**Remark:** The matrices are currently returned as a list of list of pari ``POLMOD`` objects. In the future, the ptolemy module should return the matrices as sage matrices over a `sage NumberField <http://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field.html>`_.

.. _ptolemy-example-traces:

Compute the traces
------------------

**Remark:** Requires SnapPy 2.3.2 or later.

We can compute the traces of these matrices::

    >>> sol = Manifold("m003").ptolemy_variety(2).retrieve_solutions(verbose=False)[0]
    >>> from snappy.ptolemy.matrix import matrix_trace
    >>> matrix_trace(sol.evaluate_word('a'))
    Mod(-1, x^2 - x - 1)
    >>> matrix_trace(sol.evaluate_word('b'))
    Mod(-x -1, x^2 -x -1)
    >>> matrix_trace(sol.evaluate_word('ab'))
    Mod(-x + 2, x^2 + x + 1)
    >>> matrix_trace(sol.evaluate_word('ba'))
    Mod(-x + 2, x^2 + x + 1)

**Remark:** Since this representation is irreducible, it is uniquely determined up to conjugacy by the above 4 traces, see Slide 30 of 
`Marc Culler's slides <http://www.math.illinois.edu/GEAR/resources/Culler/Culler-lecture3-slides.pdf>`_.

.. _ptolemy-examples-trace-field:

Compute the trace field for a PSL(*2*, **C**)-representation
------------------------------------------------------------

    >>> sol = Manifold("m003").ptolemy_variety(2).retrieve_solutions(verbose=False)[0]
    >>> sol.number_field()
    x^2 + x + 1

This is the Ptolemy field which is equal to the trace field if *N*\ = 2 by results of [GGZ2014]_.

.. _ptolemy-example-volume:

Compute the volume
------------------

We can also compute the volume of the representations::

    >>> sol = Manifold("m003").ptolemy_variety(2).retrieve_solutions(verbose=False)[0]
    >>> sol.volume_numerical()
    [0.E-38, 1.88266550875941 E-14]

Recall that we had an algebraic solution in the number field with defining polynomial x\ :sup:`2`\ +x+1. This number field has two embeddings into **C**, yielding two representations. This is why the result is a list of two volumes. In this case, they are both zero up to numerical precision.

.. _ptolemy-example-increase-precision:

Increase precision
------------------

We can get higher precision be setting it in pari (in decimal digits)::

    >>> sol = Manifold("m011").ptolemy_variety(2).retrieve_solutions(verbose=False)[0]
    >>> sol.volume_numerical()
    [-4.30211422042248 E-16, -0.942707362776931, 0.942707362776931]
    >>> pari.set_real_precision(40)
    15
    >>> sol.volume_numerical()
    [-1.5819817649675358086 E-40,
     -0.9427073627769277209212996030922116475902,
     0.9427073627769277209212996030922116475902]

**Remark:** This is not using interval arithmetics (although this is planned for the future). For now, the computed value of a quantity might differ from the real value by far more than the number of displayed digits suggests. To be confident about the result, we can increase the precision and see how many digits of the result are stabilizing.

.. _ptolemy-example-obstruction-class:

Ptolemy varieties for PSL(*N*, **C**)-representations
-----------------------------------------------------

The representations of ``m003`` we detected so far had trivial volume and thus cannot include the geometric representation. This is because the geometric representation is a boundary-unipotent PSL(2, **C**)-representation but not a :ref:`boundary-unipotent SL(2, C)-representation <ptolemy-boundary-unipotent>` and we only detect the latter ones above.

We can obtain the Ptolemy varieties for all :ref:`obstruction classes <obstruction-class>` to find the PSL(*N*, **C**)-representation that do not lift to boundary-unipotent SL(*N*, **C**)-representations as well::

    >>> M = Manifold("m003")
    >>> M.ptolemy_variety(N = 2, obstruction_class = 'all')
    [Ptolemy Variety for m003, N = 2, obstruction_class = 0
        c_0011_0 * c_0101_0 + c_0011_0^2 - c_0101_0^2
        c_0011_0 * c_0101_0 + c_0011_0^2 - c_0101_0^2
        - 1 + c_0011_0,
     Ptolemy Variety for m003, N = 2, obstruction_class = 1
        - c_0011_0 * c_0101_0 - c_0011_0^2 - c_0101_0^2
        - c_0011_0 * c_0101_0 - c_0011_0^2 - c_0101_0^2
        - 1 + c_0011_0]

The first Ptolemy variety in this list always corresponds to the trivial obstruction class. Let us try the non-trivial obstruction class::
    
    >>> p = M.ptolemy_variety(2, 'all')[1]
    >>> sols=p.retrieve_solutions(verbose=False)
    >>> sols.volume_numerical()
    [[2.02988321281931, -2.02988321281931]]

We now see a representation with volume twice that of a regular ideal tetrahedron. This is the geometric representation of ``m003``.
Here is python code to iterate over all obstruction classes:

    >>> for p in Manifold("m003").ptolemy_variety(2,'all'):
    ...     sols = p.retrieve_solutions(verbose=False)
    ...     print sols.volume_numerical()
    [[0.E-19, 1.88267370443418 E-14]]
    [[2.02988321281931, -2.02988321281931]]    

And in functional style::

    >>> [p.retrieve_solutions().volume_numerical() for p in Manifold("m003").ptolemy_variety(2,'all')]
    Trying to retrieve solutions from http://ptolemy.unhyperbolic.org/data/pgl2/OrientableCuspedCensus/02_tetrahedra/m003__sl2_c0.magma_out ...
    Parsing...
    Trying to retrieve solutions from http://ptolemy.unhyperbolic.org/data/pgl2/OrientableCuspedCensus/02_tetrahedra/m003__sl2_c1.magma_out ...
    Parsing...
    [[[0.E-19, 1.88267370443418 E-14]], [[2.02988321281931, -2.02988321281931]]]

**Remark**: As we see, it is not necessary to use named arguments ``N = 2`` and ``obstruction_class = 'all'`` for faster typing. However, for better readability of our code, we recommend to include the names.

A short cut for a PSL(*N*, **C**) Ptolemy variety
-------------------------------------------------

We have seen that ``M.ptolemy_variety(2, 'all')`` gives a Ptolemy variety for each obstruction class. We used ``M.ptolemy_variety(2, 'all')[3]`` to pick one, here the fourth, of those varieties. A shorter form of doing this is::

    >>> M = Manifold("m009")
    >>> M.ptolemy_variety(2, 3)
    Ptolemy Variety for m009, N = 2, obstruction_class = 3
        c_0011_0^2 + c_0101_0 * c_0101_1 + c_0101_1^2
        - c_0011_0^2 + c_0101_0^2 + c_0101_1^2
        - c_0011_0^2 - c_0101_0 * c_0101_1 - c_0101_1^2
        - 1 + c_0011_0
