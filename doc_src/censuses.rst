Census manifolds
================

Snappy comes with a large library of manifolds, which can be accessed
individually through the Manifold and Triangulation constructors but
can also be iterated through using the objects described on this
page.

SnapPy's iterators support several flexible methods for accessing
manifolds.  They can be sliced (i.e. restricted to subranges) either
by index or by volume.  Calling the iterator with keyword arguments
such as num_tets=1, betti=2 or num_cusps=3 returns an iterator which
is filtered by the specified conditions.  In addition these iterators
can determine whether they contain a given manifold.  They support
python's "A in B" syntax, and also provide an identify method which
will return a copy of the census manifold which is isometric to the
manifold passed as an argument.

..   automodule:: snappy
..   autodata:: OrientableCuspedCensus
..   autodata:: OrientableClosedCensus
..   autodata:: CensusKnots
..   autodata:: LinkExteriors
..   autodata:: HTLinkExteriors
..   autodata:: NonorientableCuspedCensus
..   autodata:: NonorientableClosedCensus

There are also:

.. toctree::
   :maxdepth: 1

   platonic_census

As instances of subclasses of ManifoldTable, the objects above
support the following methods.

..   autoclass:: snappy.database.ManifoldTable
     :members:
     :inherited-members:

Because of the large size of their datasets, the classes below
can only iterate through slices by index, and do not provide
the identification methods.

..   autoclass:: AlternatingKnotExteriors
     :members:
     :inherited-members:
..   autoclass:: NonalternatingKnotExteriors
     :members: 
     :inherited-members:
