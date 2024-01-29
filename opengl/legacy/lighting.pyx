def cyglSetStandardLighting():
    """
    Sets up our default OpenGL environment.
    """

    # Lighting intensities and location
    cdef GLfloat* ambient = [0.75, 0.75, 0.75, 1.0]
    cdef GLfloat* lightdiffuse = [0.8, 0.8, 0.8, 1.0]
    cdef GLfloat* lightspecular = [0.3, 0.3, 0.3, 1.0]
    # 2 units from the center, up and to the right
    # we should be able to control the light
    cdef GLfloat* lightposition0 = [0.3, 0.5, 3.0, 1.0]
    cdef GLfloat* lightposition1 = [0.3, -0.5, -3.0, 1.0]

    ## Set parameters that apply to all objects:
    # Remove hidden stuff
    glEnable(GL_DEPTH_TEST)
    # Allow transparency
    # glEnable(GL_ALPHA_TEST)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
    # Enable anti-aliasing of points lines and polygons
    glEnable(GL_POINT_SMOOTH)
    glEnable(GL_LINE_SMOOTH)
    #Below call is deprecated and causes odd behavior on some systems.
    #glEnable(GL_POLYGON_SMOOTH)
    # Use lights and materials to determine colors
    glEnable(GL_LIGHTING)
    # Make the Color command control ambient and diffuse material colors
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
    glEnable(GL_COLOR_MATERIAL)
    # Use interpolated shading (although colors are constant on faces)
    glShadeModel(GL_SMOOTH)
    # Define the counter-clockwise (outer) face to be the front.
    glFrontFace(GL_CCW)
    # Rasterize front and back Faces
    glDisable(GL_CULL_FACE)
    ## Set up lighting
    # Allow different properties on fronts and backs
    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0)
    # Compute specular reflections from the eye
    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, True)
    # Ambient light intensity for the entire scene
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient)
    # Enable two lights, with attenuation
    glEnable(GL_LIGHT0)
    glLightfv(GL_LIGHT0, GL_POSITION, lightposition0)
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightdiffuse)
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightspecular)
    glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION,  1.0)
    glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.1)
    glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.08)
    glDisable(GL_LIGHT1)
    glLightfv(GL_LIGHT1, GL_POSITION, lightposition1)
    glLightfv(GL_LIGHT1, GL_DIFFUSE, lightdiffuse)
    glLightfv(GL_LIGHT1, GL_SPECULAR, lightspecular)
    glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION,  1.0)
    glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.1)
    glLightf(GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 0.08)
    # Use the Model View Matrix
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
