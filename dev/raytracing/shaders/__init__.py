import os

from . import __path__ as _base_path

def _replace_compile_time_constants(shader_source, constants_dict):
    for name, value in constants_dict.items():
        shader_source = shader_source.replace(
            name, '%s' % value)
    return shader_source

_ideal_triangulation_shader_source = None

def get_ideal_triangulation_shader_source(constants_dict):

    global _ideal_triangulation_shader_source
    
    if _ideal_triangulation_shader_source is None:
        path = os.path.join(_base_path[0], 'fragment.glsl')
        _ideal_triangulation_shader_source = open(path).read()
        
    return _replace_compile_time_constants(
        _ideal_triangulation_shader_source,
        constants_dict)

_dirichlet_shader_source = None

def get_dirichlet_shader_source(constants_dict):
    
    global _dirichlet_shader_source

    if _dirichlet_shader_source is None:
        path = os.path.join(_base_path[0], 'dirichlet_fragment.glsl')
        _dirichlet_shader_source = open(path).read()

    return _replace_compile_time_constants(
        _dirichlet_shader_source,
        constants_dict)

