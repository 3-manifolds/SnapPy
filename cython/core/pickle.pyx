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
# These functions are limited to manifolds with at most 255 cusps
# such that the normal coordinates for the peripheral curves are assumed
# to be in the interval [-128, 127].  Also the Dehn filling coefficients
# must be integers in that interval.  If these conditions do not hold
# then the _to_string method can be used as an alternative. 
#
# Note that the byte sequences produced by pickle_triangulation
# are very likely to contain null bytes, so care must be taken
# in the code when converting between Python and C strings.

# Used to test if a double has an integer value:
from libc.math cimport floor

cdef IS_HUGE = <unsigned char> 1 << 7
cdef IS_BIG = <unsigned char> 1 << 6
cdef IS_COMPLETE = <unsigned char> 1 << 5

cdef pickle_triangulation(c_Triangulation *tri):
    """
    Pickle a Triangulation by constructing a TriangulationData structure
    and then serializing it.
    """
    cdef TriangulationData* tri_data
    cdef int i, flag, is_complete, num_cusps
    cdef unsigned char buf[7]
    cdef char filling[2]

    triangulation_to_data(tri, &tri_data)

    # Every pickle starts with b'pickle:'
    result = b'pickle:'

    # The first byte uses the two least signficant bits to indicate the
    # orientability type and other bits as flags.
    buf[0] = <unsigned char>tri_data.orientability
    # NOTE: tri_data.num_tetrahedra is a positive signed int.
    if tri_data.num_tetrahedra > 65535:
        flag = IS_HUGE
        buf[0] |= IS_HUGE
        buf[1] = tri_data.num_tetrahedra >> 24
        buf[2] = (tri_data.num_tetrahedra >> 16) & 0xff
        buf[3] = (tri_data.num_tetrahedra >> 8) & 0xff
        buf[4] = tri_data.num_tetrahedra & 0xff
        i = 5
    elif tri_data.num_tetrahedra > 255:
        flag = IS_BIG
        buf[0] |= IS_BIG
        buf[1] = tri_data.num_tetrahedra >> 8
        buf[2] = tri_data.num_tetrahedra & 0xff
        i = 3
    else:
        flag = 0
        buf[1] = tri_data.num_tetrahedra &0xff
        i = 2
    buf[i] = tri_data.num_or_cusps
    buf[i+1] = tri_data.num_nonor_cusps
    num_cusps = tri_data.num_or_cusps + tri_data.num_nonor_cusps
    is_complete = IS_COMPLETE
    for j in range(num_cusps):
        M, L = <double>tri_data.cusp_data[j].m, <double>tri_data.cusp_data[j].l
        if M != 0.0 or L != 0.0:
            is_complete = 0
        if M != floor(M) or L != floor(L):
            raise ValueError, 'Manifold must be pickled with _to_string.'
        if M < -128 or M > 127 or L < -128 or L > 127:
            raise ValueError, 'Manifold must be pickled with _to_string.'
    buf[0] |= is_complete
    result += buf[:i+2]
    if not is_complete:
        # Add the Dehn filling coefficients
        for j in range(num_cusps):
            # The quad double library doesn't support casting qd_real to char.
            M, L = <double>tri_data.cusp_data[j].m, <double>tri_data.cusp_data[j].l 
            filling[0] = <char>M
            filling[1] = <char>L
            result += filling[:2]
    for j in range(tri_data.num_tetrahedra):
        result += pickle_tetrahedron_data(&tri_data.tetrahedron_data[j], flag)
    return result + bytes(tri_data.name)

cdef pickle_tetrahedron_data(c_TetrahedronData* data, int flag):
    cdef int i, j, k, v, f
    cdef Py_ssize_t n = 0
    cdef unsigned short mask, bit
    cdef int x, count, curve
    cdef char curve_buf[16]
    cdef unsigned char buffer[16]

    result = bytes()
    # If this tet is in a huge or big triangulation then we use 4 or 2
    # bytes in network order to represent the indices of neighbors.
    if flag & IS_HUGE:
        for i in range(4):
            x = data.neighbor_index[i]
            buffer[n] = (x >> 24)
            buffer[n + 1] = (x >> 16) & 0xff
            buffer[n + 2] = (x >> 8) & 0xff
            buffer[n + 3] = x & 0xff
            n += 4
    elif flag & IS_BIG:
        for i in range(4):
            x = data.neighbor_index[i]
            buffer[n] = (x >> 8)
            buffer[n + 1] = x & 0xff
            n += 2
    else:
        for i in range(4):
            buffer[n] = data.neighbor_index[i]
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
            count = 0
            bit = 1
            for v in range(4):
                for f in range(4):
                    curve = data.curve[i][j][v][f]
                    if curve != 0:
                        mask |= bit
                        curve_buf[count] = <char>curve
                        count += 1
                    bit = bit << 1
            buffer[0] = mask >> 8
            buffer[1] = mask & 0xff
            # Trust Cython to copy both chars and unsigned chars into our
            # bytes object without changing any bits.
            result += buffer[:2]
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
    cdef unsigned char* p
    cdef TriangulationData tri_data
    cdef int i, num_tets, num_cusps, M, L
    cdef int n = 0

    header, tri_pickle = pickle.split(b':', 1)
    assert header == b'pickle', 'Invalid pickle byte sequence'
    p = tri_pickle
    tri_data.solution_type = not_attempted
    tri_data.volume = <Real> 0.0
    tri_data.orientability = <c_Orientability>(p[0] & 0x3)
    n = 1
    if p[0] & IS_HUGE:
        num_tets = (<int>p[n]) << 24 | (<int>p[n+1]) << 16 | (<int>p[n+2]) << 8 | (<int>p[n+3])
        n += 4
    elif p[0] & IS_BIG:
        num_tets = (<int>p[n]) << 8 | (<int>p[n+1])
        n += 2
    else:
        num_tets = p[n]
        n += 1
    tri_data.num_tetrahedra = num_tets
    tri_data.num_or_cusps = p[n]
    tri_data.num_nonor_cusps = p[n+1]
    tri_data.CS_value_is_known = 0
    tri_data.CS_value = <Real>0.0
    n += 2

    num_cusps = tri_data.num_or_cusps + tri_data.num_nonor_cusps
    if num_cusps > 0:
        # Use malloc (not mymalloc) to allocate memory for the cusp
        # data.  We free the memory before returning.
        cusps = <c_CuspData*>malloc(num_cusps*sizeof(c_CuspData))
        if cusps:
            if p[0] & IS_COMPLETE:
                for i in range(tri_data.num_or_cusps):
                    cusps[i].topology = torus_cusp
                    cusps[i].m = <Real>0.0
                    cusps[i].l = <Real>0.0
                for i in range(tri_data.num_or_cusps, num_cusps):
                    cusps[i].topology = Klein_cusp
                    cusps[i].m = <Real>0.0
                    cusps[i].l = <Real>0.0
            else:
                for i in range(tri_data.num_or_cusps):
                    M, L, n = p[n], p[n+1], n + 2
                    cusps[i].topology = torus_cusp
                    cusps[i].m = <Real>M
                    cusps[i].l = <Real>L
                for i in range(tri_data.num_or_cusps, num_cusps):
                    M, L, n = p[n], p[n+1], n + 2
                    cusps[i].topology = Klein_cusp
                    cusps[i].m = <Real>M
                    cusps[i].l = <Real>L
        else:
            raise RuntimeError, 'Failed to allocate memory for cusps'
    tri_data.cusp_data = cusps

    # Use malloc (not mymalloc) to allocate memory for the tetrahedra
    # data.  We free the memory before returning.
    tets = <c_TetrahedronData*>malloc(num_tets*sizeof(c_TetrahedronData))
    if tets == NULL:
        free(cusps)
        raise RuntimeError, 'Failed to allocate memory for tets.'
    for i in range(tri_data.num_tetrahedra):
        n = unpickle_tetrahedron_data(&tets[i], tri_pickle, n, p[0])
    tri_data.tetrahedron_data = tets

    # Create a new bytes object py_name containing the triangulation name,
    py_name = tri_pickle[n:]
    # then set tri_data.name to point to the internal buffer of py_name.
    tri_data.name = py_name
    # Hopefully, Cython takes care of the details.

    data_to_triangulation(&tri_data, &tri)
    free(tets)
    free(cusps)
    return tri

cdef int unpickle_tetrahedron_data(
    c_TetrahedronData *data, bytes pickle, int start, int flag) except *:
    cdef int i, j, v, f
    cdef int n = start
    cdef unsigned char perm
    cdef unsigned short mask, bit
    # Cython magic to allow using the internal buffer of the bytes object.
    cdef unsigned char* p = pickle

    # Get the neighbors.
    if flag & IS_HUGE:
        for i in range(4):
            data.neighbor_index[i] = ((<int>p[n]) << 24 | (<int>p[n+1]) << 16 |
                                      (<int>p[n+2]) << 8 | (<int>p[n+3]))
            n += 4
    elif flag & IS_BIG:
        for i in range(4):
            data.neighbor_index[i] = (<int>p[n] << 8) | (<int>p[n+1])
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
                        # Curve data is signed.  The effect of casting an
                        # unsigned char to a char is implementation-defined.
                        # This will fail if any bits are changed.
                        data.curve[i][j][v][f] = <char>p[n]
                        n += 1
                    else:
                        data.curve[i][j][v][f] = 0
                    bit = bit << 1

    # Return the position of the next char.
    return n
