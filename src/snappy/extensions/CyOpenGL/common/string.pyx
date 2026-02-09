def get_gl_string(string):
    cdef const char* str
    gl_string_enums = {
        'GL_VENDOR': GL_VENDOR,
        'GL_RENDERER': GL_RENDERER,
        'GL_VERSION': GL_VERSION,
        'GL_EXTENSIONS': GL_EXTENSIONS,
        'GL_SHADING_LANGUAGE_VERSION': GL_SHADING_LANGUAGE_VERSION}
    if not glGetString(GL_VERSION):
        raise RuntimeError(
            'GL strings not available - is there a current OpenGL context?')
    if string not in gl_string_enums:
        raise ValueError(
            "Invalid GL string. Must be one of %s." % ', '.join(
                [ "'%s'" % k for k in sorted(gl_string_enums.keys()) ]))
    str = <const char*>glGetString(gl_string_enums[string])
    return str.decode('ascii') if str  else 'undefined'
