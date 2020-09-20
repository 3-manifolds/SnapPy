Verified computations for closed hyperbolic 3-manifolds
=======================================================

This is an implementation of the algorithm described in the paper
of the same name (`arXiv:1904.12095 <http://arxiv.org/abs/1904.12095>`_).

It requires `SageMath <http://www.sagemath.org/>`_ (version 8.8 or later)
and `SnapPy <http://snappy.computop.org>`_ (version 2.6.1 or later).
For better performance, it also requires a particular version of
`Orb <https://github.com/DamianHeard/orb>`_ (see instructions to install Orb
below).

An example of how it can be used::

  >>> from veriClosed import *
  >>> from snappy import Triangulation
  >>> T = Triangulation('m007(3,1)').filled_triangulation()
  >>> h = compute_verified_hyperbolic_structure(T, source = 'new')
  >>> h.edge_lengths
  ...
    
This will give the edge length intervals for a verified hyperbolic
structure on the finite triangulation ``T``.

To get to the corresponding edge classes, use ``h.mcomplex.Edges`` (where
``h.mcomplex`` is an instance of ``snappy.snap.t3mlite.mcomplex.Mcomplex``).

**Notes:**

* ``filled_triangulation`` is not deterministic, thus the above code might
  fail sometimes since the finite triangulation produced by SnapPy might
  not admit a geometric structure (or such a structure could not be found).
* calling ``T.randomize()`` several times should eventually produce a
  triangulation such that we can find a geometric structure.
* Once a finite triangulation for which a geometric structure can be computed
  has been found, one can compute its isomorphism signature with
  ``T.triangulation_isosig(decorated=False)`` so that one can recreate it with
  ``Triangulation(isosig, remove_finite_vertices=False)`` in the future.
* One can obtain higher precision intervals by specifying ``bits_prec``.
  Sometimes it is necessary to increase precision for the verification of a
  hyperbolic structure to succeed.
* If ``compute_verified_hyperbolic_structure`` fails to find a hyperbolic
  structure, it raises an exception subclassed from
  ``veriClosed.verificationError.VerficationError``. Please report any
  exception raised by this function that is not of that kind: it means that
  there is an implementation error somewhere in ``veriClosed``.
* The above method uses unverified numerically methods to find an approximate
  solution to the edge equations (which we will refer to as the pre-step)
  before performing Steps I-V of the algorithm described in the paper to
  verify the hyperbolic structure. The method used in the pre-step is
  specified by the ``source`` argument.

For understanding the code in detail, the file
``veriClosed/computeVerifiedHyperbolicStructure.py`` is the best place to
start.

Pre-step method to find an unverified hyperbolic structure
----------------------------------------------------------

The following methods are available by specifying the optional argument
``source`` to ``compute_verified_hyperbolic_structure``:

* ``new``: a python-only implementation
* ``orb``: uses Orb. This works better (i.e., is faster and works for more
  triangulations so requires less randomization tries) but requires a
  particular version of Orb to be installed, see instructions below.
  (More precisely, it invokes a wrapper program
  ``orb_solution_for_snappea_finite_triangulation`` to use Orb's kernel).
* Alternatively, a path to a ``.vgm`` file can be passed as ``source`` argument.
  Such a path can be created by calling
  ``orb_solution_for_snappea_finite_triangulation`` for a ``.tri`` file
  containing a finite ``snappy.Triangulation`` saved with ``save``.



Install Orb
-----------

Some (`commit 429ca83 <https://github.com/unhyperbolic/orb/tree/429ca83>`_)
modifications are needed so that Orb can be used from ``veriClosed`` and
such that it compiles under Linux and Mac OS X. To compile this version of Orb and make it available to ``veriClosed``, do the following in some directory different from veriClosed:

* ``git clone https://github.com/unhyperbolic/orb/``
* ``cd orb``
* ``git checkout 429ca83``
* ``cd snappea/code``
* ``make``
* copy ``snappea/orb_solution_for_snappea_finite_triangulation`` into
  ``veriClosed/orb``
