To Do List
==========

- GUI

- snappy

  - Add in remaining missing features:

    - Symmetry groups
    
      - Add in actions on cusps (also to is_isometric_to)	

    - Chern-Simons 

- Documentation

  - flesh out snappy part
    
  - Installation instructions	
    
    - Linux	 
    - Windows	 
    - Sage

  - More hyperlinks
  - Make screencast

- Webpage 

  - Update webpage with t3m's fine new products, deprecate old ones.  

  - Perhaps even create an attractive modern site!

  - Idea: In public_html we will have separate directories for each
    program (SnapPy, plink, etc.) Each directory will have a
    subdirectory doc which is cloned from the documentation on a nightly
    basis.  The front page will no have anything specific about
    installation, etc.  This all will be handled by the main documenation.  

- Tasks for a later day:
   
  - Splittings 

  - Starting with generators of a lattice in PSL(2,C), building a
    fundamental domain and getting a triangulation of the corresponding
    manifold.  (Useful for building arithmetic hyperbolic 3-manifolds.)

  - dual_curves should really cache it's result and have this used by
    drill
  
  - Abelian group should he able to take any input and put it in
    canonical form, rather than simply insisting it be that way already. 
    (Cf  kernel_code/abelian_group.c/compress_abelian_group())

  - One should be able to convert a SymmetryGroup to a Sage permutation group.   

  - Also, the SymmetryGroup presentation function should be wrapped.
  There is code for this in the old SnapPeaPython.  
