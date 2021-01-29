Step-by-step examples: Part 2
=============================
    
.. _ptolemy-example-smart-lists:

The Ptolemy list type
---------------------

Recall that ``ptolemy_variety`` with ``obstruction_class='all'`` returns a list of varieties, one for each obstruction class::
    
    >>> M = Manifold("m003")
    >>> M.ptolemy_variety(2, obstruction_class = 'all')
    [Ptolemy Variety for m003, N = 2, obstruction_class = 0
        c_0011_0 * c_0101_0 + c_0011_0^2 - c_0101_0^2
        c_0011_0 * c_0101_0 + c_0011_0^2 - c_0101_0^2
        - 1 + c_0011_0,
     Ptolemy Variety for m003, N = 2, obstruction_class = 1
        - c_0011_0 * c_0101_0 - c_0011_0^2 - c_0101_0^2
        - c_0011_0 * c_0101_0 - c_0011_0^2 - c_0101_0^2
        - 1 + c_0011_0]

Also recall that ``retrieve_solutions`` was a method of a ``PtolemyVariety``. Assume we want to call ``retrieve_solutions`` for each Ptolemy variety. As in the previous example, we could write a loop such as::
 
    >>> [p.retrieve_solutions(verbose=False) for p in M.ptolemy_variety(2, 'all')]

The ptolemy module allows to do this in a much shorter way::

    >>> M.ptolemy_variety(2, 'all').retrieve_solutions(verbose=False)
    [[PtolemyCoordinates(
          {'c_0011_0': 1,
           'c_0011_1': -1,
           'c_0101_0': Mod(x, x^2 - x - 1),
	   ...,
           's_3_1': 1},
          is_numerical = False, ...)],
     [PtolemyCoordinates(
          {'c_0011_0': 1,
           'c_0011_1': -1,
           'c_0101_0': Mod(x, x^2 + x + 1),
	   ...,
           's_3_1': 1},
          is_numerical = False, ...)]]

This behavior is specific to the ptolemy module. It works with many methods of the ptolemy module that
can potentially return more than one object. These methods return a special kind of list (usually
``MethodMappingList``, a subclass of python ``list``) that tries to call the method of the given name (here ``retrieve_solutions``) with
the given arguments (here ``verbose=False``) on each element in the list (here the two Ptolemy varieties).

Since ``retrieve_solutions`` itself actually returns a list, the result is a list of lists of solutions which are of type ``PtolemyCoordinates``. The first level groups the solutions by obstruction class. The inner lists contain the different (non-Galois conjugate) solutions for each obstruction class (here, for ``m003``, each inner lists contains only one element).

Using the Ptolemy list type recursively
---------------------------------------

The list type described in the previous example works recursively. Recall that an algebraic solution to a Ptolemy variety (of type ``PtolemyCoordinates``) has a method ``volume_numerical`` that returns a list of volumes::

     >>> M=Manifold("m003")
     >>> p=M.ptolemy_variety(2,1)
     >>> sol=p.retrieve_solutions(verbose=False)[0]
     >>> sol.volume_numerical()
     [0.E-19, 1.88267370443418 E-14]

We can chain these commands together to retrieve the volumes of all boundary-unpotent PSL(2, **C**) (that are :ref:`generically decorated <ptolemy-generically-decorated>` with respect to the triangulation) in just one line::

    >>> Manifold("m003").ptolemy_variety(2,'all').retrieve_solutions(verbose=False).volume_numerical()
    [[[0.E-19, 1.88267370443418 E-14]], [[2.02988321281931, -2.02988321281931]]]

Note that the volumes of the representations are in a list of lists of lists. At the first level the volumes are grouped by obstruction class, then by Galois conjugacy.

**Remark:** There might be an extra level for witness points.

**Remark:** Unfortunately, this is not compatible with tab-autocompletion, see :ref:`later <ptolemy-example-missing-auto-completion>`.

A comparison of ``m003`` and ``m004``
-------------------------------------

We can now compare the set of volumes of ``m003`` and ``m004``:

    >>> Manifold("m003").ptolemy_variety(2,'all').retrieve_solutions(verbose=False).volume_numerical()
    [[[0.E-19, 1.88267370443418 E-14]], [[2.02988321281931, -2.02988321281931]]]
    >>> Manifold("m004").ptolemy_variety(2,'all').retrieve_solutions(verbose=False).volume_numerical()
    [[], [[-2.02988321281931, 2.02988321281931]]]

We see that the two manifolds are distinguished by their volumes of boundary-unipotent representations: ``m004`` has no representation with trivial volume (this is not a proof as in theory, there could be such a representation which is not :ref:`generically decorated <ptolemy-generically-decorated>` with respect to the given triangulation) and no representation that can be lifted to a boundary-unipotent SL(2, **C**)-representation.

A non-hyperbolic example
------------------------

We can also compute the volumes for a manifold that might be non-hyperbolic, here the complement of the 5\ :sub:`1` knot::

    >>> Manifold("5_1").ptolemy_variety(2,'all').retrieve_solutions(verbose=False).volume_numerical()
    [[], [[1.52310839130992 E-14, 0.E-37]]]

Note that one of the Ptolemy varieties is non-empty which proves that all edges of the triangulation are essential. We also see that all volumes are 0 and thus smaller than the volume 2.029883... of the figure-eight knot complement that is proven to be the smallest volume of any orientable cusped manifold. Thus, it follows from Theorem 1.3 and Remark 1.4 of [GGZ2014]_ that 5\ :sub:`1` is not hyperbolic.

**Remark:** The ptolemy module does not (yet) support interval arithmetics, otherwise, this would be a proof that 5\ :sub:`1` is not hyperbolic.


Flattening nested structures
----------------------------

If we want to loose some of the groupping, we can call ``flatten`` on the results. Here the grouping by obstruction class is lost::

    >>> Manifold("m003").ptolemy_variety(2,'all').retrieve_solutions(verbose=False).volume_numerical().flatten()
    [[0.E-19, 1.88267370443418 E-14], [2.02988321281931, -2.02988321281931]]

And now, the groupping by Galois conjugacy is lost as well, resulting in a flat list::

    >>> Manifold("m003").ptolemy_variety(2,'all').retrieve_solutions(verbose=False).volume_numerical().flatten(2)
    [0.E-19, 1.88267370443418 E-14, 2.02988321281931, -2.02988321281931]

So the result is just a flat list.

**Remark:** We cannot `overflatten`. If we give an even larger argument to flatten, the result will just stay a flat list.

.. _ptolemy-example-missing-auto-completion:

Lack of tab-autocompletion for nested structures
-------------------------------------------------

Unfortunately, the autocompletion does not list all the desired results when we have a nested structure. For example::

    >>> sols = Manifold("m003").ptolemy_variety(2,'all').retrieve_solutions(verbose=False)
    >>> sols.

When we now hit the tab key::

    >>> sols.
    sols.append   sols.extend   sols.index    sols.pop      sols.reverse  
    sols.count    sols.flatten  sols.insert   sols.remove   sols.sort

... we only get ``list`` methods, but not the desired ``volume_numerical``. One way to discover the available methods is to pick a leaf of the nested structure and hit the tab key::

    >>> sol = sols.flatten(100)[0]
    >>> sol.
    sol.N                                   sol.keys
    sol.check_against_manifold              sol.long_edge
    ...
    sol.itervalues                          sol.volume_numerical

The overview diagram might also be helpful.

Converting exact solutions into numerical solutions
---------------------------------------------------

We can turn exact solutions into numerical solutions by calling ``numerical``::

    >>> sol = Manifold("m003").ptolemy_variety(2, 1).retrieve_solutions()[0]
    >>> sol
    PtolemyCoordinates(
        {'c_0011_0': 1,
         'c_0011_1': -1,
         'c_0101_0': Mod(x, x^2 + x + 1),
	 ...
         's_3_1': 1},
        is_numerical = False, ...)
    >>> sol.numerical()
    [PtolemyCoordinates(
         {'c_0011_0': 1,
          'c_0011_1': -1,
          'c_0101_0': -0.500000000000000 - 0.866025403784439*I,
	  ...,
          's_3_1': 1},
         is_numerical = True, ...),
     PtolemyCoordinates(
         {'c_0011_0': 1,
          'c_0011_1': -1,
          'c_0101_0': -0.500000000000000 + 0.866025403784439*I,
	  ...,
          's_3_1': 1},
         is_numerical = True, ...)]

Note that the one exact (algebraic) solution turns into a list of numerical solutions which are Galois conjugates.

**Remark:** This uses the current pari precision. See the :ref:`above example <ptolemy-example-increase-precision>`, in particular, the comment about interval arithmetics.

**Remark:** Calling ``numerical()`` on a numerical solution does nothing.

**Remark:** ``CrossRatios`` also support ``numerical``.

.. _ptolemy-example-numerical-matrix:

Working with exact vs numerical solutions
-----------------------------------------

Most methods such as ``evaluate_word`` or ``cross_ratios`` work just the same way on an exact solution::

   >>> exact_sol = Manifold("m004").ptolemy_variety(2, 1).retrieve_solutions()[0]
   >>> exact_sol
   PtolemyCoordinates(
       {'c_0011_0': 1,
        'c_0011_1': -1,
        'c_0101_0': 1,
        'c_0101_1': Mod(x, x^2 + x + 1),
        ...,
	's_3_1': -1},
       is_numerical = False, ...)
   >>> exact_sol.evaluate_word('a')
   [[Mod(-2*x, x^2 + x + 1), Mod(-x - 1, x^2 + x + 1)],
    [Mod(x, x^2 + x + 1), Mod(x + 1, x^2 + x + 1)]]

... as they do on a numerical solution::

   >>> numerical_sol = sol.numerical()[0]
   >>> numerical_sol
   PtolemyCoordinates(
       {'c_0011_0': 1,
        'c_0011_1': -1,
        'c_0101_0': 1,
        'c_0101_1': -0.500000000000000 - 0.866025403784439*I,
	...,
	's_3_1': -1},
       is_numerical = False, ...)
   >>> numerical_sol.evaluate_word('a')
   [[1.00000000000000 + 1.73205080756888*I,
     -0.500000000000000 + 0.866025403784439*I],
    [-0.500000000000000 - 0.866025403784439*I,
     0.500000000000000 - 0.866025403784439*I]]

Methods with postfix ``_numerical`` are special: when applied to an exact solution, they implicitly convert it to a list
of Galois conjugate numerical solutions first. ``volume_numerical`` is an example (because volume is a transcendental function)::

    >>> exact_sol.volume_numerical()
    [-2.02988321281931, 2.02988321281931]
    >>> numerical_sol.volume_numerical()
    -2.02988321281931

.. _ptolemy-example-retrieve-numerical-solutions:

Computing numerical solutions directly
--------------------------------------

We can also directly compute numerical solutions::

    >>> M = Manifold("m004")
    >>> sols = M.ptolemy_variety(2,'all').retrieve_solutions(numerical = True)
    [[],
     [[PtolemyCoordinates(
           {'c_0011_0': 1.00000000000000 + 0.E-19*I,
            'c_0011_1': -1.00000000000000 + 0.E-19*I,
            'c_0101_0': 1.00000000000000 + 0.E-19*I,
            'c_0101_1': -0.500000000000000 - 0.866025403784439*I,
	    ...,
            's_3_1': -1},
           is_numerical = True, ...),
       PtolemyCoordinates(
           {'c_0011_0': 1.00000000000000 + 0.E-19*I,
            'c_0011_1': -1.00000000000000 + 0.E-19*I,
            'c_0101_0': 1.00000000000000 + 0.E-19*I,
            'c_0101_1': -0.500000000000000 + 0.866025403784439*I,
	    ...,
            's_3_1': -1},
           is_numerical = True, ...)]]]    

The structure is as described earlier, a list of lists of lists: first solutions are grouped by obstruction class, then by Galois conjugacy.

The advantage over going through the exact solutions is that it might be much faster
(because it can avoid computing the number field from the lexicographic Groebner basis, see later). For example, many PSL(3, **C**) examples only work when using ``numerical = True``.

.. _ptolemy-example-cross-ratios:

Computing cross ratios from Ptolemy coordinates
-----------------------------------------------

Given exact or numerical solutions to the Ptolemy variety, we can also compute the cross ratios/shape parameters::

    >>> sols = Manifold("m004").ptolemy_variety(2,'all').retrieve_solutions(verbose=False)
    >>> zs=sols.cross_ratios()
    >>> zs
    [[],
     [CrossRatios({'z_0000_0': Mod(x + 1, x^2 + x + 1),
                   'z_0000_1': Mod(x + 1, x^2 + x + 1),
                   'zp_0000_0': Mod(x + 1, x^2 + x + 1),
                   'zp_0000_1': Mod(x + 1, x^2 + x + 1),
                   'zpp_0000_0': Mod(x + 1, x^2 + x + 1),
                   'zpp_0000_1': Mod(x + 1, x^2 + x + 1)},
		  is_numerical = False, ...)]]

**Remark**: The shapes will be given as element in the Ptolemy field with defining polynomial being the second argument to ``Mod(..., ...)``, here, x\ :sup:`2`\ +x+1. The Ptolemy field is a (possibly trivial) extension of the shape field. For *N* =2, the Ptolemy field is the trace field [GGZ2014]_ and an iterated square extension of the shape field which is the invariant trace field for a cusped manifold.

And numerically, so that we can compare to SnapPy's shapes::

    >>> zs.numerical()
    [[],
     [[CrossRatios(
           {'z_0000_0': 0.500000000000000 - 0.866025403784439*I,
            'z_0000_1': 0.500000000000000 - 0.866025403784439*I,
            'zp_0000_0': 0.500000000000000 - 0.866025403784439*I,
            'zp_0000_1': 0.500000000000000 - 0.866025403784439*I,
            'zpp_0000_0': 0.500000000000000 - 0.866025403784439*I,
            'zpp_0000_1': 0.500000000000000 - 0.866025403784439*I},
           is_numerical = True, ...),
       CrossRatios(
           {'z_0000_0': 0.500000000000000 + 0.866025403784439*I,
            'z_0000_1': 0.500000000000000 + 0.866025403784439*I,
            'zp_0000_0': 0.500000000000000 + 0.866025403784439*I,
            'zp_0000_1': 0.500000000000000 + 0.866025403784439*I,
            'zpp_0000_0': 0.500000000000000 + 0.866025403784439*I,
            'zpp_0000_1': 0.500000000000000 + 0.866025403784439*I},
           is_numerical = True, ...)]]]
    >>> Manifold("m004").tetrahedra_shapes('rect')
    [0.5000000000 + 0.8660254038*I, 0.5000000000 + 0.8660254038*I]

The result is of type ``CrossRatios`` and assigns z as well as z'=1/(1-z) and z''=1-1/z a value.

.. _ptolemy-non-zero-dim-comp:

The dimension of a component
----------------------------

A Ptolemy variety might have positively dimensional components (note that this might or might not be a positively dimensional family of representations, see :ref:`here <ptolemy-generically-decorated>`). For example, the Ptolemy variety for ``m371`` and the trivial obstruction class has a 1-dimensional component. This is indicated by::

    >>> M.ptolemy_variety(2).retrieve_solutions()
    [NonZeroDimensionalComponent(dimension = 1)]

Or::

    >>> M=Manifold("m371")
    >>> M.ptolemy_variety(2).retrieve_solutions()
    [[ PtolemyCoordinates(
           {'c_0011_0': 1,
            'c_0011_1': -1,
            'c_0011_2': -1,
            'c_0011_3': Mod(-x - 1, x^2 + x + 2),
	    ...,
            's_3_4': 1},
           is_numerical = False, ...) 
       (witnesses for NonZeroDimensionalComponent(dimension = 1, free_variables = ['c_0110_2'])) ]]

The latter actually also provides a sample point (:ref:`witness <ptolemy-example-find-witness>` which we will use :ref:`later <ptolemy-example-non-zero-dim-rep>` to determine whether this corresponds to a 1-dimensional family of representations or not) on the 1-dimensional component. A ``NonZeroDimensionalComponent`` as well as ``PtolemyCoordinates`` (that correspond to 0-dimensional components of the Ptolemy variety)) has a ``dimension`` attribute, so we can do::
 
    >>> M=Manifold("m371")
    >>> sols = M.ptolemy_variety(2,'all').retrieve_solutions()
    >>> sols.dimension
    [[1], [], [0], []]

This means that the Ptolemy variety for the trivial obstruction class has a 1-dimensional component and that the Ptolemy variety of one of the other obstruction classes a 0-dimensional component.

A ``NonZeroDimensionalComponent`` is actually again a list whose elements will be witness points if witnesses have been computed for this Ptolemy variety.

**Warning:** This implies that if we ``flatten`` too much, the reported dimension becomes 0 which is the dimension of the witness point instead of 1::

    >>> sols.flatten()
    [1, 0]
    
Too much ``flatten``::
    
    >>> sols.flatten()
    [0, 0]

The advantage is that we can still call methods such as ``volume_numerical`` and actually see the volume of a witness point (it is known that the volume stays constant on a component of boundary-unipotent representations, so one witness point can tell us the volume of all representation in that component)::

    >>> sols.volume_numerical()
    [[[ [0.E-38, 0.E-38] (witnesses for NonZeroDimensionalComponent(dimension = 1, free_variables = ['c_0110_2'])) ]],
     [],
     [[4.75170196551790,
       -4.75170196551790,
       4.75170196551790,
       -4.75170196551790,
       1.17563301006556,
       -1.17563301006556,
       1.17563301006556,
       -1.17563301006556]],
     []]
