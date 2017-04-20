#  Symmetry_group

cdef class SymmetryGroup(object):
    """
    A SymmetryGroup is a group of self-isometries of hyperbolic
    3-manifold.  Instantiate as follows:

    >>> M = Manifold('m004')
    >>> M.symmetry_group()
    D4
    """
    cdef c_SymmetryGroup *c_symmetry_group
    cdef readonly _is_full_group
    cdef readonly _owns_c_symmetry_group
    
    def __cinit__(self, is_full_group, owns_c_symmetry_group):
        self.c_symmetry_group = NULL 
        self._is_full_group = is_full_group
        self._owns_c_symmetry_group = owns_c_symmetry_group

    def __dealloc__(self):
        if self._owns_c_symmetry_group:
            free_symmetry_group(self.c_symmetry_group)

    cdef _set_c_symmetry_group(self, c_SymmetryGroup * c_symmetry_group):
        if c_symmetry_group is NULL:
            raise ValueError('You tried to create an *empty* SymmetryGroup.')
        self.c_symmetry_group = c_symmetry_group

    def is_full_group(self):
        """
        Return whether the full symmetry group has been found.

        >>> S = Manifold('m004').symmetry_group()
        >>> S.is_full_group()
        True
        """
        return self._is_full_group

    def __repr__(self):
        if self.is_full_group():
            thePretext = ''
        elif self.order() == 1:
            return 'unknown'
        else:
            thePretext = 'at least '
        if self.is_abelian():
            theText = repr(self.abelian_description())
        elif self.is_dihedral():
            theText = 'D%d'%(self.order()//2)
        elif self.is_polyhedral():
            theText = self.polyhedral_description()
        elif self.is_S5():
            theText = 'S5'
        elif self.is_direct_product():
            theText =     '%s x %s' % self.direct_product_description()
        else:
            theText = 'nonabelian group of order %d'%self.order()
        
        return thePretext + theText
    
    def order(self):
        """
        Return the order of the symmetry group

        >>> S = Manifold('s000').symmetry_group()
        >>> S.order()
        4
        """
        return symmetry_group_order(self.c_symmetry_group)
    
    def is_abelian(self):
        """
        Return whether the symmetry group is abelian.

        >>> S = Manifold('m004').symmetry_group()
        >>> S.is_abelian()
        False
        """
        cdef c_AbelianGroup* abelian_description = NULL
        ans = B2B(symmetry_group_is_abelian(
                self.c_symmetry_group, &abelian_description))
        return ans

    def abelian_description(self):
        """
        If the symmetry group is abelian, return it as an AbelianGroup

        >>> S = Manifold('v3379').symmetry_group()
        >>> S.abelian_description()
        Z/2 + Z/2 + Z/2
        """
        cdef c_AbelianGroup* A
        cdef int n
        is_abelian = B2B(symmetry_group_is_abelian(self.c_symmetry_group, &A))
        if not is_abelian:
            raise ValueError('The symmetry group is not abelian.')

        coeffs = []
        for n from 0 <= n < A.num_torsion_coefficients:
                coeffs.append(A.torsion_coefficients[n])

        # Don't need to free A as it is attached to the symmetry group object
        return AbelianGroup(elementary_divisors=coeffs)
            
    
    def is_dihedral(self):
        """
        Return whether the symmetry group is dihedral.
        
        >>> S = Manifold('m004').symmetry_group()
        >>> S.is_dihedral()
        True
        """
        return B2B(symmetry_group_is_dihedral(self.c_symmetry_group))
    
    def is_polyhedral(self):
        """
        Returns whether the symmetry group is a (possibly binary)
        polyhedral group.
        """
        return B2B(symmetry_group_is_polyhedral(self.c_symmetry_group,
                                                NULL, NULL, NULL, NULL))
    
    def polyhedral_description(self):
        """
        If the symmetry group is a (possibly binary)
        polyhedral group, return a description of it.  
        """
        cdef Boolean is_binary_group
        cdef int p,q,r
        
        if not self.is_polyhedral():
            raise ValueError('The symmetry group is not polyhedral.')

        symmetry_group_is_polyhedral(self.c_symmetry_group,
                                     &is_binary_group, &p, &q, &r)

        assert p == 2

        if q == 2:
            assert(is_binary_group)
            name = 'binary dihedral group <2,2,%d>' % r
        elif q== 3:
            name = 'binary ' if is_binary_group else ''
            name += {3: 'tetrahedral group',
                     4: 'octahedral group',
                     5: 'icosohedral group'}[r]

        return name 
    
    def is_S5(self):
        """
        Returns whether the group is the symmetric group on five things.  
        """
        return B2B(symmetry_group_is_S5(self.c_symmetry_group))
    
    def is_direct_product(self):
        """
        Return whether the SymmetryGroup is a nontrivial direct
        product with at least one nonabelian factor.  
        
        >>> S = Manifold('s960').symmetry_group()
        >>> S.is_direct_product()
        True
        >>> S
        Z/4 x D3
        """
        return B2B(symmetry_group_is_direct_product(self.c_symmetry_group))
    
    def direct_product_description(self):
        """
        If the SymmetryGroup is a nontrivial direct product with at
        least one nonabelian factor, return a pair of SymmetryGroups
        consisting of the (two) factors.

        >>> S = Manifold('s960').symmetry_group()
        >>> S.direct_product_description()
        (Z/4, D3)
        """
        if not self.is_direct_product():
            raise ValueError('The symmetry group is not a nontrivial, '
                             'nonabelian direct product.')

        cdef c_SymmetryGroup* c_factor_0
        cdef c_SymmetryGroup* c_factor_1
        cdef SymmetryGroup factor_0
        cdef SymmetryGroup factor_1
        
        c_factor_0 = get_symmetry_group_factor(self.c_symmetry_group, 0)
        c_factor_1 = get_symmetry_group_factor(self.c_symmetry_group, 1)
        
        factor_0 = SymmetryGroup(True, False)
        factor_1 = SymmetryGroup(True, False)
        factor_0._set_c_symmetry_group(c_factor_0)
        factor_1._set_c_symmetry_group(c_factor_1)
        return (factor_0, factor_1)
    
    def is_amphicheiral(self):
        """
        Return whether the manifold has an orientation reversing symmetry.

        >>> S = Manifold('m004').symmetry_group()
        >>> S.is_amphicheiral()
        True
        """        
        return B2B(symmetry_group_is_amphicheiral(self.c_symmetry_group))
    
    def is_invertible_knot(self):
        """
        Return whether a one-cusped has a symmetry that acts on the
        cusp via the matrix -I.

        >>> S = Manifold('m015').symmetry_group()
        >>> S.is_invertible_knot()
        True
        """
        return B2B(symmetry_group_invertible_knot(self.c_symmetry_group))
        
    def commutator_subgroup(self):
        """
        Return the commutator subgroup of the SymmetryGroup

        >>> S = Manifold('m004').symmetry_group()
        >>> S
        D4
        >>> S.commutator_subgroup()
        Z/2
        """
        cdef c_SymmetryGroup* c_comm_subgroup
        cdef SymmetryGroup comm_subgroup

        c_comm_subgroup = get_commutator_subgroup(self.c_symmetry_group)
        comm_subgroup = SymmetryGroup(self.is_full_group(), True)
        comm_subgroup._set_c_symmetry_group(c_comm_subgroup)
        return comm_subgroup

    def abelianization(self):
        """
        Return the abelianization of the symmetry group

        >>> S = Manifold('m004').symmetry_group()
        >>> S.abelianization()
        Z/2 + Z/2
        """
        
        if not self.is_full_group():
            raise ValueError('The full symmetry group is not known.')

        cdef c_SymmetryGroup* c_abelianization
        cdef SymmetryGroup abelianization

        c_abelianization = get_abelianization(self.c_symmetry_group)
        abelianization = SymmetryGroup(self.is_full_group(), True)
        abelianization._set_c_symmetry_group(c_abelianization)
        return abelianization.abelian_description()

    def center(self):
        """
        Return the center of the symmetry group

        >>> S = Manifold('m004').symmetry_group()
        >>> S.center()
        Z/2
        """
        
        if not self.is_full_group():
            raise ValueError('The full symmetry group not known.')

        cdef c_SymmetryGroup* c_center
        cdef SymmetryGroup center

        c_center = get_center(self.c_symmetry_group)
        center = SymmetryGroup(self.is_full_group(), True)
        center._set_c_symmetry_group(c_center)
        return center.abelian_description()

    def multiply_elements(self, i, j):
        """
        Returns the product of group elements i and j.  The convention
        is that products of symmetries read right to left.  That is,
        the composition (symmetry[i] o symmetry[j]) acts by first
        doing symmetry[j], then symmetry[i].

        >>> S = Manifold('m004').symmetry_group()
        >>> S.multiply_elements(2, 3)
        1
        """
        cdef int prod
        order = self.order()
        for x in [i,j]:
            if not (0 <= x < order):
                raise ValueError('The symmetry group has only %d '
                                 'elements.' % order)

        return symmetry_group_product(self.c_symmetry_group, i, j)

    def isometries(self):
        """
        Return a detailed list of all the isometries in the symmetry group.

        >>> S = Manifold('s959').symmetry_group()
        >>> isoms = S.isometries()
        >>> isoms[8]
        0 -> 1   1 -> 0 
        [-1 -1]  [ 0  1]
        [ 1  0]  [-1 -1]
        Does not extend to link
        """
        cdef IsometryList *isometries

        isometries = get_symmetry_list(self.c_symmetry_group)
        ans = IsometryListToIsometries(isometries)
        return ans
