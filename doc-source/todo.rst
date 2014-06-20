To Do List
==========

- GUI
  
  - Allow you to get the info (position, radius) of a horoball by clicking it.  

- Documentation

  - More hyperlinks

  - Expand tutorial 

  - Make more screencasts

- snappy 
  
  - Longitude detection in manifolds with a single torus cusp (ie, the slope that dies in homology).  Requested by Saul.  

- Kernel 

  - Merge Jeff's changes into SnapPy's copy of the kernel.

  - After this, consider making SnapPy's copy the canonical one by
    definition. 
 
  - Improve the situation with uFatalErrors.  

- Minor 

  - dual_curves should really cache it's result and have this used by
    drill
  
  - One should be able to convert a SymmetryGroup to a Sage permutation group.   

  - Also, the SymmetryGroup presentation function should be wrapped.
    There is code for this in the old SnapPeaPython.  

- Ambitious

  - A new basic (sub)class: S3Knot (and/or S3Link).
 
  - Consider adding a HeegaardSplitting class 

  - Consider merging our t3m project and normal surface code into
    SnapPy. 

  - Redo much of Snap in the context of Sage/SnapPy.   

     - Add a method for computing tetrahedron shapes to arbitrary precision.

     - Add methods for computing invariant trace fields and related number
       fields.

     - Add a method which implements and extends Harriet Moser's
       algorithm to allow SnapPy to prove that a manifold is hyperbolic.

