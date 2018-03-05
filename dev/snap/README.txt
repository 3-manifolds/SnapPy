> I remember you saying that snappy.snap had trouble with some >
> examples that Snap could do, but I saw no sign of that.

That was a problem in the past, but Matthias recently made some
changes to the algebraic number recognition code so this may well no
longer be the case. In theory, snappy.snap should actually be better
since it's using a more sophisticated LLL implementation (fpLLL vs
PARI). I guess one test would be whether we can easily replicate the
computations that are in "snap-pari/snap_data", in particular the
contents of "closed.fields" and "cusped.fields".

> What do you think we would need to do before claiming that SnapPy
> can do anything Snap can do?

This is not a complete list, but:

(a) Computation of invariant quaternion algebra and determining its structure (which primes ramify etc).

(b) Determine if all traces are algebraic integers.

(c) Determining if given manifold is arithmetic (easy from (a) and (b), I think).

(d) Snap typically gets the holonomy representation into a smaller field than SnapPy, since they do some clever initial conjugation.

(e) Easily access the cusp field.

(f) Compute the borel regulator, i.e. the list of volumes of all the Galois conjugates of the holonomy representation.

I don't think any of this is super-difficult, especially given that we
have the Snap source. (Sage's support of quaternion algebras over
fields other than Q is not complete, but it does have some of what
would be needed already.)

Things not directly related to number theory that Snap can do but
SnapPy can't. These are just ordinary precision C++ code, and so are
things that we could probably import into SnapPy without too much
difficulty if we wanted.

(i) Compute words in fundamental group associated to elements of the
length spectrum

(ii) Compute ortholines, i.e. distances between lifts of the same
closed geodesic. Other related functions related to tubes about
geodesics.

(iii) Find hidden symmetries

(iv) Compute the commensurator 

(v) Compute eta invariant

(vi) Alterative algorithm for finding symmetry groups which provides
additional information (e.g. matrices realizing the symmetries).

Best,

Nathan

