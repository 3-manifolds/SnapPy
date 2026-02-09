import os

from . import __path__ as _base_path

def get_texture_paths():
    return [ os.path.join(_base_path[0], name)
             for name in ['NonGeometric.png',
                          'Eye.png'] ]

def _replace_compile_time_constants(shader_source, constants_dict):
    for name, value in constants_dict.items():
        shader_source = shader_source.replace(
            name, ('%s' % value).encode())
    return shader_source


_triangulation_shader_source = None


def get_triangulation_shader_source_and_ubo_descriptors(
        constants_dict,
        defs_dict = {}):

    global _triangulation_shader_source

    if _triangulation_shader_source is None:
        path = os.path.join(_base_path[0], 'fragment.glsl')
        _triangulation_shader_source = open(path, 'rb').read()

    src = _replace_compile_time_constants(
        _triangulation_shader_source,
        constants_dict)

    if defs_dict:
        header, footer = src.split(b'\n', 1)
        def_block = (
            '\n\n'
            '// { Compile time defines\n')
        for k, v in sorted(defs_dict.items()):
            def_block += '#define %s' % k
            if not v is None:
                def_block += ' %r' % v
            def_block += '\n'
        def_block += '// } Compile time defines\n\n'
        src = header + def_block.encode() + footer

    num_geodesic_segments = defs_dict.get('num_geodesic_segments', 0)
    num_additional_horospheres = defs_dict.get('num_additional_horospheres', 0)
    num_eyeballs = defs_dict.get('num_eyeballs', 0)
    num_tets = constants_dict[b'##num_tets##']

    uniform_block_names_sizes_and_offsets = [
        ('TetrahedraCombinatorics',
         (64 + 64) * num_tets,
         { 'otherTetNums' : 0,
           'otherFaceNums' : 64 * num_tets }),
        ('TetrahedraBasics',
         (64 + 64 + 256) * num_tets,
         { 'R13Vertices' : 0,
           'planes' : 64 * num_tets,
           'SO13tsfms' : (64 + 64) * num_tets } ),
        ('Colors',
         (64 + 96 + 64) * num_tets,
         { 'face_color_indices' : 0,
           'edge_color_indices' : 64 * num_tets,
           'vertex_color_indices' : (64 + 96) * num_tets } ),
        ('TetrahedraEdges',
         192 * num_tets,
         { 'R13EdgeEnds' : 0 } ),
        ('TetCuspMatrices',
         (256 + 256) * num_tets,
         { 'tetToCuspMatrices' : 0,
           'cuspToTetMatrices': 256 * num_tets } ),
        ('MargulisTubes',
         (64 + 64) * num_tets,
         { 'margulisTubeTails': 0,
           'margulisTubeHeads' : 64 * num_tets}),
        ('geodesics',
         (16 + 16 + 16 + 16) * num_geodesic_segments + 16 * (num_tets + 1),
         { 'geodesicTails': 0,
           'geodesicHeads': 16 * num_geodesic_segments,
           'geodesicIndex': (16 + 16) * num_geodesic_segments,
           'geodesicTubeRadiusParam': (16 + 16 + 16) * num_geodesic_segments,
           'geodesicOffsets': (16 + 16 + 16 + 16) * num_geodesic_segments }),
        ('eyeballs',
         (16 + 64 + 64) * num_eyeballs + 16 * (num_tets + 1),
         { 'eyeballPositions' : 0,
           'eyeballInvEmbeddings' : 16 * num_eyeballs,
           'eyeballEmbeddings' : (16 + 64) * num_eyeballs,
           'eyeballOffsets' : (16 + 64 + 64) * num_eyeballs }),
        ('additionalHorospheres',
         (16 + 16) * num_additional_horospheres + 16 * (num_tets + 1),
         { 'horosphereVec': 0,
           'horosphereCuspIndex': 16 * num_additional_horospheres,
           'horosphereOffsets' : (16 + 16) * num_additional_horospheres}),
        ('edgeMidpoints',
         16 * 6 * num_tets,
         { 'edgeMidpointVec': 0 })
         ]

    return src, uniform_block_names_sizes_and_offsets
