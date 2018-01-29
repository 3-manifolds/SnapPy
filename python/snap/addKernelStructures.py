"""
This file contains methods to attach data from SnapPea kernel to a
t3mlite.Mcomplex.
"""

from .t3mlite import simplex

Infinity = "Infinity"

def _clean_ideal_vertices(vertices):
    """
    The SnapPea kernel gives us a large number for infinity.
    Convert it to infinity.
    """
    return [ x if abs(x) < 10**20 else Infinity for x in vertices]

def addShapes(mcomplex, shapes):
    """
    Given a shape for each tetrahedron, add to the tetrahedron.

    The only assumption made here about the type of shapes is that
    they support the basic arithmetic operations. In particular,
    they can be SnapPy numbers or complex intervals.
    """
    
    for tet, z in zip(mcomplex.Tetrahedra, shapes):
        zp  = 1 / (1 - z)
        zpp = (z - 1) / z
        tet.ShapeParameters = {
            simplex.E01: z,
            simplex.E23: z,
            simplex.E02: zp,
            simplex.E13: zp,
            simplex.E03: zpp,
            simplex.E12: zpp
            }
        
def addChooseGeneratorInfo(mcomplex, choose_generators_info):
    """
    Expects a Mcomplex and the result of Manifold._choose_generator_info().

    Adds GeneratorsInfo to each tetrahedron. This encodes generator_index and
    generator_status of a SnapPea Triangulation as described in
    choose_generators.c. However, we only store one number, its absolute value
    giving the generator_index and its sign the generator_status.

    We also set ChooseGenInitialTet of the Mcomplex to be what SnapPea
    considers the base tetrahedron when computing the vertices of the
    fundamental domain.

    We also add the vertices of a fundamental domain as given by the SnapPea
    kernel as SnapPeaIdealVertices to each tetrahedron. We care about these
    numbers when orienting the base tetrahedron to have consistency with the
    SnapPea kernel.
    """
    
    for tet, info in zip(mcomplex.Tetrahedra, choose_generators_info):
        tet.SnapPeaIdealVertices = dict(
            zip(simplex.ZeroSubsimplices,
                _clean_ideal_vertices(info['corners'])))
        tet.GeneratorsInfo = dict(
            zip(simplex.TwoSubsimplices,
                info['generators']))
        if info['generator_path'] == -1:
            mcomplex.ChooseGenInitialTet = tet

def reindexCuspsAndAddPeripheralCurves(
    mcomplex, cusp_indices_and_peripheral_curve_data):
    """
    Expects a Mcomplex and the result of
    Manifold._get_cusp_indices_and_peripheral_curve_data().

    It rearranges the Vertices of the mcomplex to match the ordering
    of the cusps in the SnapPea kernel and adds the peripheral curves
    in a format analogous to the kernel.
    """

    cusp_indices, curves = cusp_indices_and_peripheral_curve_data

    def process_row(curves):
        return {
            vertex : {
                face : curves[4 * i + j]
                for j, face in enumerate(simplex.TwoSubsimplices) }
            for i, vertex in enumerate(simplex.ZeroSubsimplices) }
        
    for i, tet in enumerate(mcomplex.Tetrahedra):

        tet.PeripheralCurves = [
            # meridian
            [ process_row(curves[4 * i + 0]),   # right-handed sheet
              process_row(curves[4 * i + 1]) ], # left-handed sheet
            # longitude
            [ process_row(curves[4 * i + 2]),   # right-handed sheet
              process_row(curves[4 * i + 3]) ]] # left-handed sheet
        
        for vertex, cusp_index in zip(simplex.ZeroSubsimplices,
                                      cusp_indices[i]):
            tet.Class[vertex].Index = cusp_index
            
    mcomplex.Vertices.sort(key = lambda vertex : vertex.Index)

    for index, vertex in enumerate(mcomplex.Vertices):
        if not index == vertex.Index:
            raise Exception("Inconsistencies with vertex indices")
