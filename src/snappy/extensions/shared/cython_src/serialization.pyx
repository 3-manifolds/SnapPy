toml_header = """\
class = "{class_name}"
name = "{name}"

[computed_data]
solution_type = "{solution_type}"
"""

cdef get_class_name(obj):
    return str(obj.__class__).split("'")[1].split(".")[-1]

cdef solution_type_values = {
    not_attempted: 'not_attempted',
    geometric_solution: 'geometric_solution',
    nongeometric_solution: 'nongeometric_solution',
    flat_solution: 'flat_solution',
    degenerate_solution: 'degenerate_solution',
    other_solution: 'other_solution',
    no_solution: 'no_solution',
    externally_computed: 'externally_computed',
    orb_partially_flat_solution: 'orb_partially_flat_solution',
    orb_step_failed: 'orb_partially_flat_solution',
    orb_invalid_solution: 'orb_invalid_solution',
}
    
cusp_topology_values = [
    'sphere',
    'projective plane',
    'torus',
    'Klein bottle',
    'orientable higher genus',
    'nonorientable higher  genus',
    'unknown']

toml_cusp = """
[[cusp]]
topology = \"{topology}\"
filling = {filling}
"""

cdef get_cusp(c_CuspData *data):
    return toml_cusp.format(
        topology=cusp_topology_values[<int>data.topology],
        # The Cython compiler crashes without these casts to double.
        filling=str([<double>data.m, <double>data.l]))

toml_tetrahedron = """
[[tetrahedron]]
neighbors = {neighbors}
gluings = {gluings}
cusps = {cusps}
{curves}"""

cdef get_tetrahedron(c_TetrahedronData *tet):
    cdef c_TetrahedronData data = tet[0]
    neighbors = [data.neighbor_index[i] for i in range(4)]
    gluings = [[data.gluing[i][j] for j in range(4)] for i in range(4)]
    cusps = [data.cusp_index[i] for i in range(4)]
    meridian_right = [[data.curve[0][0][i][j] for i in range(4)] for j in range(4)]
    meridian_left = [[data.curve[0][1][i][j] for i in range(4)] for j in range(4)]
    longitude_right = [[data.curve[1][0][i][j] for i in range(4)] for j in range(4)]
    longitude_left = [[data.curve[1][1][i][j] for i in range(4)] for j in range(4)]
    curves = ''
    if any(sum(meridian_right, [])):
        curves += 'curves.meridian.right = ' + str(meridian_right) + '\n'
    if any(sum(meridian_left, [])):
        curves += 'curves.meridian.left = ' + str(meridian_left) + '\n'
    if any(sum(longitude_right, [])):
        curves += 'curves.longitude.right = ' + str(longitude_right) + '\n'
    if any(sum(longitude_left, [])):
        curves += 'curves.longitude.left = ' + str(longitude_left) + '\n'
    return toml_tetrahedron.format(
        neighbors=str(neighbors),
        gluings=str(gluings),
        cusps=str(cusps),
        curves=curves)

toml_shape = """\
[[computed_data.shape]]
real = {real}
imag = {imag}
"""

cdef snappy_object_to_string(Triangulation obj):
    """Serialize a SnapPy object in Toml format.
    The class of the object must inherit from Triangulation
    """
    cdef c_Triangulation *c_triangulation = obj.c_triangulation
    cdef TriangulationData *data
    triangulation_to_data(c_triangulation, &data)
    cdef double volume = <double> data.volume
    cdef double cs = <double> data.CS_value
    result = ''
    # write the header
    class_name = get_class_name(obj)
    name = data.name.decode('utf8')
    solution_type = solution_type_values[data.solution_type]
    result += toml_header.format(class_name=class_name,
                                 name=name,
                                 solution_type=solution_type)
    if (data.solution_type == geometric_solution or
        data.solution_type == nongeometric_solution):
        result += f'volume = {volume}\n'
        for i in range(data.num_tetrahedra):
            real = <double>data.tetrahedron_data[i].filled_shape.real
            imag = <double>data.tetrahedron_data[i].filled_shape.imag
            result += toml_shape.format(real=real, imag=imag)
    if data.CS_value_is_known:
        result += f'charn_simons = {cs}\n';
    # write the cusps
    for i in range(data.num_or_cusps + data.num_nonor_cusps):
        result += get_cusp(&data.cusp_data[i])
    # write the tetrahedra
    for i in range(data.num_tetrahedra):
        result += get_tetrahedron(&data.tetrahedron_data[i])
    return result

def save_snappy_object(Triangulation obj, str filename):
    """Write a Serialized SnapPy object to a Toml file.
    The class of the object must inherit from Triangulation
    """
    with open(filename, 'w') as toml_file:
        toml_file.write(snappy_object_to_string(obj))
