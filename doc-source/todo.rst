To Do List
==========

- GUI

- snappy

  - Switch from DirichletDomain(M) to M.dirichlet_domain()
  - Same for CuspNeighborhood(M), even though it is not finished.
  - Add doctests interface update for FundamentalGroup, etc.
  - Same for the censuses.  
  - Make sure things still work under Sage 4.0
  - Add in remaining missing features:

    - Symmetry groups (including action on cusps)
    - Setting the peripheral curves
    - Chern-Simons 

- Documentation

  - flesh out snappy part
    
    - Should methods of Manifolds inherited unchanged from
      Triangulation appear also in the latter's page?

  - Installation instructions	
    
    - Linux	 
    - Windows	 

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
  
