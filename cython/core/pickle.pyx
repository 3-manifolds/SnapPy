# Cython functions for pickling and unpickling Triangulations
#
# The functions pickle_triangulation and unpickle triangulation
# perform fairly straightforward serialization and deserialization of
# the kernel's TriangulationData structure, via Python bytes objects.
#
# For complete orientable manifolds, this process preserves all
# combinatorial properties of the triangulation, including the
# orderings of the tetrahedra and of the vertices in each tetrahedron
# as well as the normal curve description of the peripheral curves.
# (The peripheral curves are not preserved for manifolds which are not
# orientable or not complete.)
#
# These functions are limited to manifolds with at most 2^15 - 1
# tetrahedra and at most 255 cusps.  Also, the normal coordinates for
# the peripheral curves are assumed to be in the interval [-128, 127]
#
# Note that the byte sequences produced by pickle_triangulation
# are very likely to contain null bytes, so care must be taken
# in the code when converting between Python and C strings.

cdef pickle_triangulation(c_Triangulation *tri):
    """
    Pickle a Triangulation by constructing a TriangulationData structure
    and then serializing it.
    """
    cdef TriangulationData* tri_data
    cdef int i, is_big, use_cusp_data
    cdef char buf[5]

    triangulation_to_data(tri, &tri_data)

    # Every pickle starts with b'pickle:'
    result = b'pickle:'

    is_big = 1 if tri_data.num_tetrahedra > 127 else 0
    # The first byte indicates the orientability type and whether the triangulation
    # is big.
    buf[0] = <char>tri_data.orientability
    if is_big:
        # A big triangulation is flagged by setting bit 7 of the first byte
        # number as the first byte of the pickle.
        buf[0] |= 0x80
        buf[1] = <char>(tri_data.num_tetrahedra >> 8)
        buf[2] = <char>(tri_data.num_tetrahedra & 0xff)
        i = 3
    else:
        buf[1] = <char>tri_data.num_tetrahedra
        i = 2

    # We only use the cusp data if all cusps are complete.  Otherwise, when
    # unpickling we set both cusp counts to 0 so the kernel will build its own
    # cusps.
    use_cusp_data = True
    for j in range(tri_data.num_or_cusps):
            if (tri_data.cusp_data[j].m != 0.0 or
                tri_data.cusp_data[j].l != 0.0):
                use_cusp_data = False
                break
    if use_cusp_data:
        buf[i] = tri_data.num_or_cusps
        buf[i+1] = tri_data.num_nonor_cusps
    else:
        buf[i] = buf[i+1] = 0
    result += buf[:i+2]
    for i in range(tri_data.num_tetrahedra):
        result += pickle_tetrahedron_data(&tri_data.tetrahedron_data[i],
                                          is_big)
    return result + bytes(tri_data.name)

cdef pickle_tetrahedron_data(c_TetrahedronData* data, int is_big):
    cdef int i, j, k, v, f
    cdef Py_ssize_t n = 0
    cdef short x, mask, bit
    cdef int count, curve
    cdef char curve_buf[16]
    cdef char buffer[16]

    result = bytes()
    # If this tet is in a big triangulation then we use 2 bytes in
    # network order to represent the indices of its neighbors.
    if is_big:
        for i in range(4):
            x = <short>data.neighbor_index[i]
            buffer[n] = <char>(x >> 8)
            buffer[n + 1] = <char>(x & 0xff)
            n += 2
    else:
        for i in range(4):
            buffer[n] = <char>data.neighbor_index[i]
            n += 1

    # Pack each permutation into a byte, like SnapPea does.
    cdef int *perm
    for i in range(4):
        perm = &data.gluing[i][0]
        buffer[n] = perm[0] | perm[1]<<2 | perm[2]<<4 | perm[3]<<6
        n += 1

    # Use 1 byte per cusp index.  This assumes < 256 cusps.
    for i in range(4):
        buffer[n] = <unsigned char>data.cusp_index[i]
        n += 1
    result += buffer[:n]

    # The curve data consists of 64 integers, which are often mostly
    # zeros.  We optimize a bit to take advantage of that feature.  We
    # use chars for the integers.  We divide the data into 4 groups of
    # 16 and represent each group with a 16 bit mask, as 2 bytes in
    # network order, followed by one byte for each nonzero bit in the
    # mask.
    for i in range(2):
        for j in range(2):
            mask = 0
            count = 2
            bit = 1
            for v in range(4):
                for f in range(4):
                    curve = data.curve[i][j][v][f]
                    if curve != 0:
                        mask |= bit
                        curve_buf[count] = <char>curve
                        count += 1
                    bit = bit << 1
            curve_buf[0] = <char>(mask >> 8)
            curve_buf[1] = <char>(mask & 0xff)
            result += curve_buf[:count]
    return result

cdef c_Triangulation* unpickle_triangulation(bytes pickle) except *:
    """
    Unpickle a Triangulation by deserializing a TriangulationData structure
    and then converting it to a Triangulation.
    """
    cdef c_TetrahedronData* tets = NULL
    cdef c_CuspData* cusps = NULL
    cdef c_Triangulation *tri
    cdef char* p
    cdef TriangulationData tri_data
    cdef int i, num, is_big
    cdef int n = 0

    header, tri_pickle = pickle.split(b':', 1)
    assert header == b'pickle', 'Invalid pickle byte sequence'
    p = tri_pickle
    is_big = 1 if p[0] & 0x80 else 0
    tri_data.solution_type = not_attempted
    tri_data.volume = <Real> 0.0
    tri_data.orientability = <c_Orientability>(p[0] & 0x7f)
    n = 1
    if is_big:
        num = p[n] << 8 | p[n+1]
        n += 2
    else:
        num = p[n]
        n += 1
    tri_data.num_tetrahedra = num
    tri_data.num_or_cusps = p[n]
    tri_data.num_nonor_cusps = p[n+1]
    tri_data.CS_value_is_known = 0
    tri_data.CS_value = <Real>0.0
    n += 2

    cdef num_cusps = tri_data.num_or_cusps + tri_data.num_nonor_cusps
    if num_cusps > 0:
        # Use malloc (not mymalloc) to allocate memory for the cusp
        # data.  We free the memory before returning.
        cusps = <c_CuspData*>malloc(num_cusps*sizeof(c_CuspData))
        for i in range(tri_data.num_or_cusps):
            cusps[i].topology = torus_cusp
            cusps[i].m = <Real>0.0
            cusps[i].l = <Real>0.0
        for i in range(tri_data.num_or_cusps, num_cusps):
            cusps[i].topology = Klein_cusp
            cusps[i].m = <Real>0.0
            cusps[i].l = <Real>0.0
    tri_data.cusp_data = cusps

    # Use malloc (not mymalloc) to allocate memory for the tetrahedra
    # data.  We free the memory before returning.
    tets = <c_TetrahedronData*>malloc(num*sizeof(c_TetrahedronData))
    for i in range(tri_data.num_tetrahedra):
        n = unpickle_tetrahedron_data(&tets[i], tri_pickle, n, is_big)
    tri_data.tetrahedron_data = tets

    # Create a new bytes object py_name containing the triangulation name,
    py_name = tri_pickle[n:]
    # then set tri_data.name to point to the internal buffer of py_name.
    tri_data.name = py_name
    # Cython takes care of the details.

    data_to_triangulation(&tri_data, &tri)
    if tets:
        free(tets)
    if cusps:
        free(cusps)
    return tri

cdef int unpickle_tetrahedron_data(
    c_TetrahedronData *data, bytes pickle, int start, int is_big) except *:
    cdef int i, j, v, f
    cdef int n = start
    cdef unsigned char perm
    cdef unsigned short mask, bit
    # Cython magic to allow using the internal buffer of the bytes object.
    cdef char* p = pickle

    # Get the neighbors.
    if is_big:
        for i in range(4):
            data.neighbor_index[i] = <int>p[n] << 8 | <int>p[n+1]
            n += 2
    else:
        for i in range(4):
            data.neighbor_index[i] = <int>p[n]
            n += 1

    # Get the gluing permutations.
    for i in range(4):
        perm = <unsigned char>p[n]
        for j in range(4):
            data.gluing[i][j] = <int>(perm & 0x3)
            perm = perm >> 2
        n += 1

    # Get the cusp indices.
    for i in range(4):
        data.cusp_index[i] = <int>p[n]
        n += 1

    # Get the curve data
    for i in range(2):
        for j in range(2):
            bit = 1
            mask = <unsigned char>p[n] << 8 | <unsigned char>p[n+1]
            n += 2
            for v in range(4):
                for f in range(4):
                    if mask & bit:
                        data.curve[i][j][v][f] = p[n]
                        n += 1
                    else:
                        data.curve[i][j][v][f] = 0
                    bit = bit << 1

    # Return the position of the next char.
    return n
