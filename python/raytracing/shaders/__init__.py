import os

from . import __path__ as _base_path

def _replace_compile_time_constants(shader_source, constants_dict):
    for name, value in constants_dict.items():
        shader_source = shader_source.replace(
            name, ('%s' % value).encode())
    return shader_source

_ideal_triangulation_shader_source = None

def get_ideal_triangulation_shader_source_and_ubo_descriptors(constants_dict):

    global _ideal_triangulation_shader_source
    
    if _ideal_triangulation_shader_source is None:
        path = os.path.join(_base_path[0], 'fragment.glsl')
        _ideal_triangulation_shader_source = open(path, 'rb').read()

    src = _replace_compile_time_constants(
        _ideal_triangulation_shader_source,
        constants_dict)

    num_tets = constants_dict[b'##num_tets##']

    uniform_block_names_sizes_and_offsets = [
        ('TetrahedraBasics',
         (64 + 64 + 256) * num_tets,
         { 'R13Vertices' : 0,
           'planes' : 64 * num_tets,
           'SO13tsfms' : (64 + 64) * num_tets } ),
        ('TetCuspMatrices',
         (256 + 256) * num_tets,
         { 'tetToCuspMatrices' : 0,
           'cuspToTetMatrices': 256 * num_tets } ),
        ('MargulisTubes',
         (64 + 64) * num_tets,
         { 'margulisTubeTails': 0,
           'margulisTubeHeads' : 64 * num_tets}) ]
    
    return src, uniform_block_names_sizes_and_offsets

_dirichlet_shader_source = None

def get_dirichlet_shader_source(constants_dict):
    
    global _dirichlet_shader_source

    if _dirichlet_shader_source is None:
        path = os.path.join(_base_path[0], 'dirichlet_fragment.glsl')
        _dirichlet_shader_source = open(path).read()

    return _replace_compile_time_constants(
        _dirichlet_shader_source,
        constants_dict)

