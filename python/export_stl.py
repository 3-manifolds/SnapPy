
import math

def facet_stl(triangle):
    vertex1, vertex2, vertex3 = triangle
    a = (vertex3[0]-vertex1[0], vertex3[1]-vertex1[1], vertex3[2]-vertex1[2])
    b = (vertex2[0]-vertex1[0], vertex2[1]-vertex1[1], vertex2[2]-vertex1[2])
    normal = (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])
    return ''.join([
        '  facet normal %f %f %f\n' % tuple(normal),
        '    outer loop\n',
        '      vertex %f %f %f\n' % tuple(vertex1),
        '      vertex %f %f %f\n' % tuple(vertex2),
        '      vertex %f %f %f\n' % tuple(vertex3),
        '    endloop\n',
        '  endfacet\n'
        ])

def subdivide_triangles(triangles, num_subdivisions):
    if num_subdivisions == 0:
        for triangle in triangles:
            yield triangle
    elif num_subdivisions == 1:
        # A small function for getting the midpoint of two points.
        midpoint = lambda P, Q: ((P[0] + Q[0]) / 2, (P[1] + Q[1]) / 2, (P[2] + Q[2]) / 2)
        for (x, y, z) in triangles:
            yield (x, midpoint(x, y), midpoint(x, z))
            yield (midpoint(y, x), y, midpoint(y, z))
            yield (midpoint(z, x), midpoint(z, y), z)
            yield (midpoint(x, y), midpoint(y, z), midpoint(z, x))
    else:
        for triangle in subdivide_triangles(subdivide_triangles(triangles, 1), num_subdivisions - 1):
            yield triangle
    return

def projection(triangle, cutoff_radius):
    ''' Return the projection of a point in the Klein model to the Poincare model. '''
    x, y, z = triangle
    scale = min(1 / (1 + math.sqrt(max(0, 1 - (x**2 + y**2 + z**2)))), cutoff_radius)
    return (scale*x, scale*y, scale*z)


def klein_stl(face_dicts):
    ''' Yield triangles describing these faces. '''
    for face in face_dicts:
        vertices = face['vertices']
        for i in range(len(vertices)-2):
            yield (vertices[0], vertices[i+1], vertices[i+2])
    return

def klein_cutout_stl(face_dicts, shrink_factor=0.9):
    ''' Yield triangles describing these faces after removing a fraction of the interior.
    
    The fraction removed is given by shrink_factor. '''
    for face in face_dicts:
        vertices = face['vertices']
        center = [sum(vertex[i] for vertex in vertices) / len(vertices) for i in range(3)]
        new_vertices = [[vertex[i] + (center[i] - vertex[i]) / 3 for i in range(3)] for vertex in vertices]
        new_inside_points = [[point[i] * shrink_factor for i in range(3)] for point in new_vertices]
        for i in range(len(new_vertices)):
            yield (new_vertices[i], new_inside_points[(i+1) % len(new_vertices)], new_inside_points[i])
            yield (new_vertices[i], new_vertices[(i+1) % len(new_vertices)], new_inside_points[(i+1) % len(new_vertices)])
            yield (vertices[i], new_vertices[(i+1) % len(vertices)], new_vertices[i])
            yield (vertices[i], vertices[(i+1) % len(vertices)], new_vertices[(i+1) % len(vertices)])
            # We have to go in the opposite direction this time as the normal should point in towards O.
            yield tuple(tuple(shrink_factor * coord for coord in point) for point in (vertices[i], new_vertices[i], new_vertices[(i+1) % len(vertices)]))
            yield tuple(tuple(shrink_factor * coord for coord in point) for point in (vertices[i], new_vertices[(i+1) % len(vertices)], vertices[(i+1) % len(vertices)]))
    return

def poincare_stl(face_dicts, num_subdivisions=5, cutoff_radius=0.9):
    ''' Yield the output of klein_stl(face_dicts, ...) after applying projection to every vertex produced. '''
    for triangle in subdivide_triangles(klein_stl(face_dicts), num_subdivisions):
        yield (projection(triangle[0], cutoff_radius), projection(triangle[1], cutoff_radius), projection(triangle[2], cutoff_radius))
    return

def poincare_cutout_stl(face_dicts, num_subdivisions=3, shrink_factor=0.9, cutoff_radius=0.9):
    ''' Yield the output of klein_cutout_stl(face_dicts, ...) after applying projection to every vertex produced. '''
    for triangle in subdivide_triangles(klein_cutout_stl(face_dicts, shrink_factor), num_subdivisions):
        yield (projection(triangle[0], cutoff_radius), projection(triangle[1], cutoff_radius), projection(triangle[2], cutoff_radius))
    return

def stl(face_dicts, model='klein', cutout=False, num_subdivisions=3, shrink_factor=0.9, cutoff_radius=0.9):
    """
    Yield the lines of an stl file corresponding to the solid given by face_dicts that is suitable for 3d printing.
    
    Arguments can be given to modify the model produced:
        model='klein' - (alt. 'poincare') the model of HH^3 to use.
        cutout=False - remove the interior of each face
        shrink_factor=0.9 - the fraction to cut out of each face
            cuttoff_radius=0.9 - maximum rescaling for projection into Poincare model
        num_subdivision=3 - number of times to subdivide for the Poincare model
    For printing domains in the Poincare model, cutoff_radius is critical for avoiding infinitely
    thin cusps, which cannot be printed.
    """
    if shrink_factor < 0 or shrink_factor > 1:
        raise ValueError('shrink_factor must be between 0 and 1.')
    
    if model == 'klein' and cutout:
        output = klein_cutout_stl(face_dicts, shrink_factor=shrink_factor)
    elif model == 'klein' and not cutout:
        output = klein_stl(face_dicts)
    elif model == 'poincare' and cutout:
        output = poincare_cutout_stl(face_dicts, num_subdivisions=num_subdivisions, shrink_factor=shrink_factor, cutoff_radius=cutoff_radius)
    elif model == 'poincare' and not cutout:
        output = poincare_stl(face_dicts, num_subdivisions=num_subdivisions, cutoff_radius=cutoff_radius)
    else:
        raise ValueError('Unknown model. Known models: \'klein\' and \'poincare\'.')
    
    yield 'solid\n'
    for triangle in output:
        yield facet_stl(triangle)
    yield 'endsolid\n'
    return

