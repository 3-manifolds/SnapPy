Internals of verified computations
==================================



Naming
------

The names of methods containing ``check`` will raise an exception if
the desired property cannot be certified. There are different types of
Exceptions to indicate how the certification failed. This type can be
used by other methods to perform some action such as changing the
triangulation or increasing precision or to give up.

The user-facing methods have names starting with ``verify`` or
``verified`` and will fail more gracefully returning ``False`` or
``None`` in such a case.



Generating certified shape intervals
------------------------------------

The recommended way to obtain certified intervals for the shapes is via
``manifold.tetrahedra_shapes(intervals=True)`` as :doc:`described
earlier <verify>`. Here we document the ``KrawczykShapesEngine`` and
``IntervalNewtonShapesEngine`` which is implemented internally to
generate the intervals. It is of interest for those users who want to
understand the underlying interval math and experiment with the Newton
interval method or the Krawczyk test.  ``CertifiedShapesEngine`` is an
alias of either ``KrawczykShapesEngine`` or
``IntervalNewtonShapesEngine`` to determine the default method used by
verify.

..   automodule:: snappy.verify
..   autoclass:: CertifiedShapesEngine
     :members:
     :inherited-members:

..   autoclass:: IntervalNewtonShapesEngine
     :members:
     :inherited-members:

..   autoclass:: KrawczykShapesEngine
     :members:
     :inherited-members:


Verification of hyperbolicity
-----------------------------

Methods containing ``check`` will raise an exception if the desired property
cannot be certified. Methods containing ``verify`` or ``verified`` will fail
more gracefully returning ``False`` or ``None`` in such a case.

..   autofunction:: snappy.verify.verifyHyperbolicity.check_logarithmic_gluing_equations_and_positively_oriented_tets

Cusp cross sections
-------------------

..   autoclass:: snappy.verify.RealCuspCrossSection
     :members:
     :inherited-members:

..   autoclass:: snappy.verify.ComplexCuspCrossSection
     :members:
     :inherited-members:

Verified canonical cell decompositions
--------------------------------------

..   autofunction:: snappy.verify.verifyCanonical.interval_checked_canonical_triangulation
..   autofunction:: snappy.verify.verifyCanonical.exactly_checked_canonical_retriangulation

Exact computations for cusp cross sections
------------------------------------------

..   automodule:: snappy.verify.squareExtensions

..   autofunction:: snappy.verify.squareExtensions.find_shapes_as_complex_sqrt_lin_combinations
..   autoclass:: SqrtLinCombination
     :members:
..   autoclass:: ComplexSqrtLinCombination
     :members:

Exceptions
----------

..   automodule:: snappy.verify.exceptions
     :members:
