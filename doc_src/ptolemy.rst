The ptolemy module
==================

..   automodule:: snappy.ptolemy
		  
What is the ptolemy module?
---------------------------

This module provides tools to find boundary-unipotent representations of an
oriented
3-manifold into PSL(*N*, **C**). It was used, for example, to generate the tables
of representations at
`ptolemy.unhyperbolic.org <http://ptolemy.unhyperbolic.org/html/summary.html>`_.

The ptolemy module can use `magma <http://magma.maths.usyd.edu.au/magma/>`_ for
the computations necessary to find the representations or it can automatically
retrieve the necessary computations from a database we provide that contains the
computations for all manifolds and *N* at
`ptolemy.unhyperbolic.org <http://ptolemy.unhyperbolic.org/html/summary.html>`_.
In particular, the database includes the computations for all orientable census
manifolds for PSL(2, **C**) and the examples below work without magma for these
manifolds.
This makes it useful for people who do not have magma (see :ref:`here <ptolemy-example-retrieve-exact-solutions>`).

The ptolemy module is still under development. Please report bugs or email suggestions to: enischte at gmail dot com.

Examples of what the Ptolemy module can do:

* Give the :ref:`equations for the reduced Ptolemy varieties<ptolemy-example-basic>` to find all :ref:`generically decorated <ptolemy-generically-decorated>` and :ref:`boundary-unipotent <ptolemy-boundary-unipotent>` representations into SL(*N*, **C**).
* Same for representations into PSL(*N*, **C**) using :ref:`Obstruction classes <ptolemy-example-obstruction-class>`.
* Retrieve the :ref:`exact <ptolemy-example-retrieve-exact-solutions>` and :ref:`numerical <ptolemy-example-retrieve-numerical-solutions>` solutions for these varieties.
* Compute the matrices for these representations :ref:`exactly <ptolemy-example-matrices>` or :ref:`numerically <ptolemy-example-numerical-matrix>` given a generator of or word in the fundamental group.
* Compute the :ref:`traces <ptolemy-example-traces>` for these matrices.
* Compute the :ref:`trace field <ptolemy-examples-trace-field>` of the representation.
* Compute the :ref:`boundary holonomoy <ptolemy-example-boundary-holonomy>`.
* Compute the corresponding :ref:`shape parameters/cross ratios <ptolemy-example-cross-ratios>` for a representation.
* Compute the :ref:`volume <ptolemy-example-volume>` and :ref:`complex volume <ptolemy-example-complex-volume>` (volume + i Chern-Simons) for a representation.
* Compute solutions to the Ptolemy varieties :ref:`using magma or sage <ptolemy-example-using-magma-sage>`.
* Find :ref:`positively dimensional components <ptolemy-non-zero-dim-comp>` of the Ptolemy variety.
* Find :ref:`positively dimensional families <ptolemy-example-non-zero-dim-rep>` of boundary-unipotent representations.

The ptolemy module tries to have powerful and flexible functionality but at the same time be easy to use. For example, we can get the volumes of boundary-unipotent PSL(2, **C**)-representations of ``m011`` in just one line::

    >>> Manifold("m011").ptolemy_variety(2,'all').retrieve_solutions().volume_numerical()
    [[[-4.30211422042248 E-16, -0.942707362776931, 0.942707362776931]],
     [[4.64255370258293 E-15, 3.94215909915729 E-15, -2.78183391239608, 2.78183391239608]]]

Here, we see that ``m011`` has a representation that is not Galois-conjugate to the geometric representation and that has the volume 0.9427... of the Weeks manifold.

**Remark:** Please check your Internet connection if the above example does not work.

Documentation
-------------

.. toctree::
   :maxdepth: 2

   ptolemy_prelim
   ptolemy_examples1
   ptolemy_examples2
   ptolemy_examples3
   ptolemy_examples4
   ptolemy_classes



