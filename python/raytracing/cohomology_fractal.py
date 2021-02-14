from snappy.ptolemy.homology import homology_basis_representatives_with_orders

"""
Helpers to compute weights for faces corresponding to a 2-cocycle. Used to
visualize cohomology fractals.

The weights for the faces are stored in the order
[face 0 tet 0, face 1 tet 0, ... face 3 tet 0, face 0, tet 1, ...].

Note that we also use face classes (consisting each of two
faces). Some helpers convert weights stored per face class to weights
stored per face per tetrahedron.
"""

def face_var_name_to_index(var_name):
    """
    Convert variable name to index in array of weights.

        >>> face_var_name_to_index('s_2_5')
        22

    """

    name, face_index, tet_index = var_name.split('_')
    if name != 's':
        raise AssertionError(
            "Variable name '%s' for face class invalid" % var_name)
    return 4 * int(tet_index) + int(face_index)

def value_for_face_class(weights, face_class):
    """
    Given weights per face per tetrahedron and a face_class as encoded
    by the ptolemy module, extract the weight for that face_class and
    perform consistency check.
    """
    
    sgn, power, repr0, repr1 = face_class

    # Weight for one representative of face class
    val0 = weights[face_var_name_to_index(repr0)]
    # Weight for other representative of face class
    val1 = weights[face_var_name_to_index(repr1)]
    
    # Check weights match
    if abs(val0 - sgn * val1) > 1e-6:
        raise ValueError("Weights for identified faces do not match")

    return val0

def check_weights_valid(trig, weights):
    """
    Given a SnapPy triangulation and weights per face per tetrahedron,
    check they are consistent across glued faces and form a 2-cocycle.

        >>> from snappy import Manifold
        >>> check_weights_valid(Manifold("m015"), [1, 0, 0, 0, -1, 0, 0, 1, -1, 0, 0, 0])
        >>> check_weights_valid(Manifold("m004"), [0, 1, 0, 1, 0, 0, -1, -2])
        Traceback (most recent call last):
           ...
        ValueError: Weights for identified faces do not match
        >>> check_weights_valid(Manifold("m003"), [1, 0, 0, 2, -1, 0, 0, -2])
        Traceback (most recent call last):
           ...
        ValueError: Weights are not a 2-cocycle.

    """

    face_classes = trig._ptolemy_equations_identified_face_classes()

    check_face_class_weights_valid(
        trig,
        [ value_for_face_class(weights, face_class)
          for face_class in face_classes ])

def check_face_class_weights_valid(trig, weights):
    """
    Given a SnapPy triangulation and weights per face class, check they
    form a 2-cocycle.
    """
    matrix, edge_labels, face_labels = (
        trig._ptolemy_equations_boundary_map_2())
    for row in matrix:
        total = sum(entry * weight for entry, weight in zip(row, weights))
        if abs(total) > 1e-6:
            raise ValueError("Weights are not a 2-cocycle.")

def compute_signs_and_face_class_indices(trig):
    """
    An array helpful to convert weights per face class to weights per
    face per tetrahedron. The entries are per face per tetrahedron and
    are a pair (sign adjustment, face class index).
    """

    result = (4 * trig.num_tetrahedra()) * [ None ]
    face_classes = trig._ptolemy_equations_identified_face_classes()
    for i, face_class in enumerate(face_classes):
        sgn, power, repr0, repr1 = face_class
        result[face_var_name_to_index(repr0)] = ( +1, i)
        result[face_var_name_to_index(repr1)] = (sgn, i)
    return result

def rational_cohomology_basis(trig):
    """
    Given a SnapPy triangulation, give 2-cocycles encoded as weights
    per face per tetrahedron that generate the second rational cohomology
    group.

        >>> from snappy import Manifold
        >>> rational_cohomology_basis(Manifold("m125"))
        [[-1, -4, -2, 0, 1, 0, 0, 0, 0, 4, 1, 2, -1, 0, 0, 0], [0, 5, 2, 0, 0, 0, 0, 1, 0, -5, 0, -2, 0, 0, 0, -1]]

    """

    signs_and_face_class_indices = compute_signs_and_face_class_indices(
        trig)

    matrix2, edge_labels, face_labels = (
        trig._ptolemy_equations_boundary_map_2())
    matrix3, face_labels, tet_labels = (
        trig._ptolemy_equations_boundary_map_3())

    two_cocycles_and_orders = (
        homology_basis_representatives_with_orders(
            matrix2, matrix3, 0))

    rational_two_cocycles = [
        two_cocycle
        for two_cocycle, order in two_cocycles_and_orders
        if order == 0 ]

    # Expand two cocycles which assign a weight to each face class
    # to weights per face per tetrahedron
    return [
        [ sgn * two_cocycle[index]
          for sgn, index in signs_and_face_class_indices ]
        for two_cocycle in rational_two_cocycles ]

def compute_weights_basis_class(trig, cohomology_class):
    """
    Convenience method to quickly access cohomology classes for
    M.inside_view().

        >>> from snappy import Manifold
        >>> compute_weights_basis_class(Manifold("m004"), None)
        (None, None, None)
        >>> compute_weights_basis_class(Manifold("m003"), [1, 0, 0, 1, -1, 0, 0, -1])
        ([1, 0, 0, 1, -1, 0, 0, -1], None, None)
        >>> compute_weights_basis_class(Manifold("m003"), 0)
        (None, [[1, 0, 0, 1, -1, 0, 0, -1]], [1.0])
        >>> compute_weights_basis_class(Manifold("m125"), 0)
        (None, [[-1, -4, -2, 0, 1, 0, 0, 0, 0, 4, 1, 2, -1, 0, 0, 0], [0, 5, 2, 0, 0, 0, 0, 1, 0, -5, 0, -2, 0, 0, 0, -1]], [1.0, 0.0])
        >>> compute_weights_basis_class(Manifold("m125"), 1)
        (None, [[-1, -4, -2, 0, 1, 0, 0, 0, 0, 4, 1, 2, -1, 0, 0, 0], [0, 5, 2, 0, 0, 0, 0, 1, 0, -5, 0, -2, 0, 0, 0, -1]], [0.0, 1.0])
        >>> compute_weights_basis_class(Manifold("m125"), [0.5, 0.5])
        (None, [[-1, -4, -2, 0, 1, 0, 0, 0, 0, 4, 1, 2, -1, 0, 0, 0], [0, 5, 2, 0, 0, 0, 0, 1, 0, -5, 0, -2, 0, 0, 0, -1]], [0.5, 0.5])

    """

    # If no cohomology_class specified, just return None
    if cohomology_class is None:
        return None, None, None

    try:
        as_list = list(cohomology_class)
    except TypeError:
        as_list = None

    if as_list:
        # User can give a weight for each tetrahedron
        if len(as_list) == 4 * trig.num_tetrahedra():
            check_weights_valid(trig, as_list)
            return as_list, None, None

        # User can specify super position of basis two-cocycles,
        # that is a number for each element in the basis.
        basis = rational_cohomology_basis(trig)
        if len(as_list) == len(basis):
            return None, basis, as_list
        
        raise ValueError(
            ("Expected array of length %d or %d either assigning one number "
             "for each basis vector of the second rational cohomology group "
             "or one weight per face "
             "and tetrahedron.") % (len(basis), 4 * trig.num_tetrahedra()))
    else:
        # User can just specify an integer to pick one of the basis
        # two-cocycles generating the cohomology.
        basis = rational_cohomology_basis(trig)

        c = range(len(basis))[cohomology_class]

        return None, basis, [ 1.0 if i == c else 0.0
                              for i in range(len(basis)) ]
