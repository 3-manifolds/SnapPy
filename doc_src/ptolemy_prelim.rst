Mathematical preliminaries
==========================

Given a triangulation, the ptolemy module will produce a system of equation that is equivalent to
the reduced Ptolemy variety (see [GTZ2011]_, Section 5 of [GGZ2012]_, and Proposition 4.7 of [GGZ2014]_).

A solution
to this system of equations where no Ptolemy coordinate is zero yields a :ref:`boundary-unipotent <ptolemy-boundary-unipotent>`
SL(*N*, **C**)-representation, respectively, PSL(*N*, **C**)-representation (see :ref:`obstruction-class`).

Note that a solution where some Ptolemy coordinates are zero might not have enough information
to recover the representation - thus the ptolemy module discards those and thus might miss some
boundary-unipotent representations for the chosen triangulation (see :ref:`ptolemy-generically-decorated`).
This is the same problem that the
gluing equations for finding PGL(2, **C**)-representations suffer from where simplices in the developing
map can be degenerate and yielding cross ratios that are 0, 1, or :math:`\infty`\ .

.. _ptolemy-boundary-unipotent:

Boundary-unipotent
------------------

We call a SL(*N*,  **C**)-representation *boundary-unipotent* if each peripheral subgroup is taken to
a conjugate of the unipotent group *P* of upper unit-triangular matrices. Similarly, we call
a PSL(*N*, **C**)-representation *boundary-unipotent* if each peripheral subrgroup is taken to a conjugate
of the unipotent group of PSL(*N*, **C**), i.e., the image of *P* under the projection SL(*N*, **C**)\ :math:`\rightarrow`\ PSL(*N*, **C**).

Note that even when boundary-unipotent PSL(*N*, **C**)-representation can be lifted to an
SL(*N*, **C**)-representation, the lift might no longer be boundary-unipotent, i.e., there might be
a peripheral curve whose image now is conjugate to an upper triangular matrix with an *N*-th root
of unity on the diagonal. For example, if the manifold is hyperbolic and has one cusp,
any lift of the geometric representation will take the longitude
to a matrix with trace -2 and is thus not boundary-unipotent as SL(2, **C**)-representation.

.. _obstruction-class:

Obstruction class
-----------------

Given a boundary-unipotent PSL(*N*, **C**)-representation, we obtain an *obstruction class* in H\ :sup:`2`\ (M,\ :math:`\partial`\ M; **Z**/*N*)
that is trivial if and only if the representations lifts to a boundary-unipotent SL(*N*, **C**)-representation (see Section 9.1 of [GTZ2011]_ and Section 1.3 of [GGZ2014]_).
Given a triangulation and an element in H\ :sup:`2`\ (M,\ :math:`\partial`\ M; **Z**/*N*), the Ptolemy variety can be modified to find
the boundary-unipotent 
PSL(*N*, **C**)-representations with that obstruction class (for *N* > 2 this is implemented here but has not been published yet).

Note that two elements in H\ :sup:`2`\ (M,\ :math:`\partial`\ M; **Z**/*N*)
related by multiplication by an element in (**Z**/*N*)\ :sup:`*` yield Ptolemy
varieties corresponding to picking different Galois conjugates for the *N*-th root of unity. Thus, it is enough
to consider a representative for each element in the quotient H\ :sup:`2`\ (M,\ :math:`\partial`\ M; **Z**/*N*)/(**Z**/*N*)\ :sup:`*`\ .

.. _ptolemy-psl-multiplicity:

SL(*N*, **C**) vs PSL(*N*, **C**)
---------------------------------

The reduced Ptolemy variety for the trivial obstruction class can have several points (say *d*) giving different SL(*N*, **C**)-representations that are the same as PSL(*N*, **C**)-representations. Similarly, for the non-trivial obstruction class we can have *d* points in the reduced Ptolemy variety yielding the same PSL(*N*, **C**)-representation.

The degree *d* of this correspondence is the size of H\ :sup:`1`\ (\ :math:`\hat{M}`\ ; **Z**/*N*) where :math:`\hat{M}` is the cell complex obtained from collapsing each cusp to a point.

.. _ptolemy-generically-decorated:

Generically decorated representations
-------------------------------------

We want to point out two important facts when using the reduced Ptolemy variety to find boundary-unipotent representations:

* We miss representations that are not generically decorated (as mentioned above). This happens but rarely.
* A positively dimensional component in the reduced Ptolemy variety might mean two things (which we can distinguish by looking at the images of the peripheral groups):
     *  a positively dimensional family of boundary-unipotent representations or
     *  a family of decorations of the same representation.

The reason for this is that the reduced Ptolmey variety does not parametrize representations but generically decorated representations (which can also be thought of as development maps). We just list the facts about decorations important to us here and refer the reader for details to Section 4 of [GTZ2011]_ and Section 8 of [GGZ2012]_ (where decoration would be called *B*-decoration or (SL(*N*), **C**), *B*)-decoration with *B* the Borel group of upper triangular matrices):

* Every boundary-unipotent representation of a cusped manifold admits a decoration. The set of decorations of a representation is intrinsic to the representation and independent of the triangulation.
* The representation determines the decoration uniquely if and only if the representation is boundary-non-degenerate (which most representations are).
* Given a decorated representation and an ideal triangulation of a cusped manifold, we obtain Ptolemy coordinates.
* If all the resulting Ptolemy coordinates are non-zero, we call the representation *generically decorated* - a notion that depends on the chosen triangulation.
* The reduced Ptolemy variety parametrizes generically decorated and boundary-unipotent representations.

  
.. _ptolemy-reduced-variety:

Reduced Ptolemy variety
-----------------------

We will actually always use the reduced Ptolemy variety, i.e., the system of equation that consists of the Ptolemy relations (always of the form
:math:`\pm` c\ :sub:`...` c\ :sub:`...` :math:`\pm` c\ :sub:`...` c\ :sub:`...`  :math:`\pm` c\ :sub:`...` c\ :sub:`...`\) and extra equations fixing an appropriate set of (N-1) Ptolemy coordinates per cusp as described in Proposition 4.7 of [GGZ2014]_. This is because the Ptolemy relations alone admit an action by (**C**\ :sup:`*`\ )\ :sup:`(N-1)` for each cusp that does not change the representation it yields.

In other words, the Ptolemy variety parametrizes *P*-decorations and the reduced Ptolemy variety parametrizes *B*-decorations.

Future work
-----------

In unpublished work, we developed an algorithm that takes some triangulation of a manifold and constructs a set of triangulations and corresponding Ptolemy varieties (with extra edge relations) such that we can guarantee that all boundary-unipotent PSL(2,C)-representations are found - not just the ones that are generically decorated with respect to the chosen triangulation. This is inspired by [S2009]_. The ptolemy module might support this in the future.

In [Z2014]_, the Ptolemy variety was extended to detect non boundary-unipotent representations as well. The ptolemy module might produce these varities in the future. This might offer another way of computing A-polynomials - that when combined with the above algorithm is guarenteed to be the full A-polynomial and not just a factor of it.

References
----------

.. [S2009] Henry Segerman: A generalisation of the deformation variety, http://arxiv.org/abs/0904.1893
.. [GTZ2011] Stavros Garoufalidis, Dylan P. Thurston, and Christian K. Zickert: The Complex Volume of SL(n,C)-Representations of 3-Manifolds, http://arxiv.org/abs/1111.2828
.. [GGZ2012] Stavros Garoufalidis, Matthias Goerner, and Christian K. Zickert: Gluing Equations for PGL(n,C)-Representations of 3-Manifolds, http://arxiv.org/abs/1207.6711
.. [GGZ2014] Stavros Garoufalidis, Matthias Goerner, and Christian K. Zickert: The Ptolemy Field of 3-Manifold Representations, http://arxiv.org/abs/1401.5542
.. [Z2014] Christian K. Zickert: Ptolemy coordinates, Dehn invariant, and the A-polynomial, http://arxiv.org/abs/1405.0025
