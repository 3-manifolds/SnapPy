import snappy

def preserves_orientation(iso):
    return all(A[0,0]*A[1,1] - A[0,1]*A[1,0] == 1 for A in iso.cusp_maps())

def are_orient_pres_isometric(M, N):
    return any(preserves_orientation(iso)
               for iso in M.is_isometric_to(N, True))

for M in snappy.LinkExteriors[1.0:]:
    if M.solution_type() != 'contains degenerate tetrahedra':
        E = snappy.Link(M.name()).exterior()
        if not are_orient_pres_isometric(M, E):
            print(M.name())
    
    
