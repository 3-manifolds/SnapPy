include "CyOpenGL.pxi"
include "CyOpenGLU.pxi"

# This is part of the UCS2 hack.  Of course, for this to work this
# file should not contain any unicode strings.

cdef public UCS2_hack (char *string, Py_ssize_t length, char *errors) :   
    raise RuntimeError, """
Don't use unicode strings in SnapPy Cython files!"""

import Tkinter as Tk_
import os, sys, platform, tkMessageBox
from colorsys import hls_to_rgb
from math import sqrt, ceil, floor
from random import random

def glVersion():
    cdef char *gl_version
    gl_version = <char *> glGetString(GL_VERSION)
    return gl_version

cdef class vector3:
    """
    A simple real 3-dimensional vector which supports addition,
    subtraction and right multiplication or division by scalars.
    Attributes include its norm and the square of its norm.
    """
    cdef readonly double x, y, z, norm_squared, norm

    def __cinit__(self, triple):
        self.x, self.y, self.z = map(float, triple)
        self.norm_squared = self.x*self.x + self.y*self.y + self.z*self.z 
        self.norm = sqrt(self.norm_squared)

    def __repr__(self):
        return '< %s, %s, %s >'%(self.x, self.y, self.z)

    def __add__(self, vector):
        return vector3([self.x+vector.x, self.y+vector.y, self.z+vector.z])

    def __sub__(self, vector):
        return vector3([self.x-vector.x, self.y-vector.y, self.z-vector.z])

    def __mul__(self, scalar):
        return vector3([self.x*scalar, self.y*scalar, self.z*scalar])

    def __div__(self, scalar):
        return vector3([self.x/scalar, self.y/scalar, self.z/scalar])

cdef class GL_context:
    """
    Sets up our default OpenGL environment.
    """

    def __cinit__(self):
        # Lighting intensities and location
        cdef float* ambient = [0.5, 0.5, 0.5, 1.0]
        cdef float* lightdiffuse = [0.8, 0.8, 0.8, 1.0]
        cdef float* lightspecular = [0.3, 0.3, 0.3, 1.0]
        # 2 units from the center, up and to the right
        # we should be able to control the light
        cdef float* lightposition = [0.1, 0.1, 1.5, 1.0]

        ## Set parameters that apply to all objects:
        # Remove hidden stuff
        glEnable(GL_DEPTH_TEST)
        # Allow transparency
        # glEnable(GL_ALPHA_TEST)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        # Enable anti-aliasing of lines and points
        glEnable(GL_POINT_SMOOTH)
        glEnable(GL_LINE_SMOOTH)
        # Use lights and materials to determine colors
        glEnable(GL_LIGHTING)
        # Make the Color command control ambient and diffuse material colors
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
        glEnable(GL_COLOR_MATERIAL)
        # Use interpolated shading (although colors are constant on faces)
        glShadeModel(GL_SMOOTH)
        # Define the counter-clockwise (outer) face to be the front.
        glFrontFace(GL_CCW);
        # Rasterize front and back Faces
        glDisable(GL_CULL_FACE);
        ## Set up lighting
        # Allow different properties on fronts and backs
        glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0)
        # Compute specular reflections from the eye
        glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, True)
        # Ambient light intensity for the entire scene
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient)
        # Enable one light, with attenuation
        glEnable(GL_LIGHT0)
        glLightfv(GL_LIGHT0, GL_POSITION, lightposition)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, lightdiffuse)
        glLightfv(GL_LIGHT0, GL_SPECULAR, lightspecular)
        glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION,  1.0)
        glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.2)
        glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.08)
        # Use the Model View Matrix
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity()

cdef class GLU_context:
    """
    Holds pointers to whatever structures are required by GLU.
    """
    cdef GLUquadric* glu_quadric

    def __cinit__(self):
        self.glu_quadric = gluNewQuadric()

    def __dealloc__(self):
        gluDeleteQuadric(self.glu_quadric)


cdef class GLobject:
    """
    Base class for the objects in our OpenGL scene.  Think of a
    GLobject as a helper who knows how to draw a certain type of
    geometrical object with specific color and material characteristics.
    The geometrical characteristics, e.g. radius or vertex locations,
    should be passed as arguments to the object's draw method.
    """
    cdef GLfloat color[4]
    cdef GLfloat front_specular[4]
    cdef GLfloat back_specular[4]
    cdef GLfloat front_shininess
    cdef GLfloat back_shininess

    def __cinit__(self, *args,
                  color = [0.8, 0.8, 0.8, 1.0],
                  front_specular = [0.8, 0.8, 0.8, 1.0], 
                  back_specular = [0.8, 0.8, 0.8, 1.0],
                  front_shininess = 0.0,
                  back_shininess = 0.0,
                  **kwargs):
        cdef int n
        for n from 0 <= n < 4:
            try:
                self.color[n] = color[n]
            except:
                print self.color[n]
            self.front_specular[n] = float(front_specular[n])
            self.back_specular[n] = float(back_specular[n])
        self.front_shininess = float(front_shininess)
        self.back_shininess = float(back_shininess)

    def set_material(self):
        glMaterialfv(GL_FRONT, GL_SPECULAR, self.front_specular)
        glMaterialf(GL_FRONT, GL_SHININESS, self.front_shininess)
        glMaterialfv(GL_BACK,  GL_SPECULAR, self.back_specular)
        glMaterialf(GL_BACK,  GL_SHININESS, self.back_shininess)
        glColor4fv(self.color)

    def draw(self, *args, **kwargs):
        """
        Issue the OpenGL commands to draw this object.
        (Override in subclasses)
        """

    def build_display_list(self, list_id, *args):
        """
        Generate a display list containing the commands to draw this object.
        (Override in subclasses)
        """

cdef class Sphere(GLobject):
    """
    Draw a sphere.  Use a wire frame when filled=False, solid otherwise.
    The sphere is drawn as a GLU quadric.
    """
    cdef GLUquadric* glu_quadric

    def __cinit__(self, *args, GLU_context GLU, **kwargs):
        self.glu_quadric = GLU.glu_quadric

    def __init__(self,
                 GLU_context GLU,
                 filled=False,
                 color=[0.8,0.8,0.8,0.3],
                 front_specular = [0.8, 0.8, 0.8, 1.0], 
                 back_specular = [0.8, 0.8, 0.8, 1.0],
                 front_shininess = 50.0,
                 back_shininess = 0.0
                 ):
        if not filled:
            gluQuadricDrawStyle(self.glu_quadric, GLU_LINE)
        else:
            gluQuadricDrawStyle(self.glu_quadric, GLU_FILL)
        gluQuadricNormals(self.glu_quadric, GLU_SMOOTH)
     
    def draw(self, GLdouble radius, GLint slices, GLint stacks):
        self.set_material()
        # We put the north pole on the y-axis. 
        glPushMatrix()
        glLoadIdentity()
        glRotatef(90, 1.0, 0.0, 0.0)
        gluSphere(self.glu_quadric, radius, slices, stacks)
        glPopMatrix()

    def build_display_list(self, list_id, radius, slices, stacks):
        glNewList(list_id, GL_COMPILE) 
        self.draw(radius, slices, stacks)
        glEndList()

class TriangleMesh:
    """
    A triangle which can tessellate itself.
    """
    def __init__(self, vertices):
        self.vertices = vertices
        self.triangles = [(0,1,2)]

    def __repr__(self):
        return str(self.triangles)

    def __getitem__(self, n):
        x, y, z = self.triangles[n]
        return (self.vertices[x], self.vertices[y], self.vertices[z])

    def subdivide(self):
        """
        Replace each triangle by four triangles:
                       z
                     /   \
                    zx - yz
                   /  \ /  \
                  x -  xy - y
        New midpoint vertices are appended to the vertex list.
        """
        new_triangles = []
        V = self.vertices
        for triangle in self.triangles:
            x, y, z = triangle
            n = len(V)
            self.vertices.append((V[x] + V[y])/2)
            self.vertices.append((V[y] + V[z])/2)
            self.vertices.append((V[z] + V[x])/2)
            xy, yz, zx = n, n+1, n+2 
            new_triangles += [(x, xy, zx), (xy, yz, zx),
                              (zx, yz, z), (xy, y, yz)]
        self.triangles = new_triangles

cdef class PoincareTriangle(GLobject):
    """
    Draws a geodesic triangle in the Poincare model.  The geometric
    parameters are the vertex coordinates in the Klein model plus the
    coordinates of the center of the sphere which represents the plane
    of the triangle in the Poincare model.  The Poincare vertices are
    constructed by projecting the Klein vertices onto the sphere from
    the center.
    The triangle is drawn as a mesh, using vertex and index arrays.
    """
    cdef vertices, center, mesh, count
    cdef GLfloat* nv_array
    cdef GLushort* indices
    cdef GLuint buffers[2]
#   cdef int useVBO
#   The switch useVBO is currently set to False by default.  The VBO
#   code worked fine with NVidia drivers in OS X but there were segfaults
#   inside the i915 drivers with linux/X11.  Performance seems fine
#   either way.  It doesn't seem worthwhile to sniff graphics cards.
#   XX In fact, the MinGW MESA library doesn't even define the VBO calls,
#   so I have now commented them out until support improves.

    def __init__(self, vertices, center, subdivision_depth=4,
#                 useVBO=False, 
                 **kwargs):
        self.vertices = vertices
        self.center = center
#        self.useVBO = useVBO
        self.mesh = TriangleMesh(vertices)
        for n in range(subdivision_depth):
            self.mesh.subdivide()
        self.build_arrays()

    def __dealloc__(self):
        free(self.nv_array)
        free(self.indices)
#        if self.useVBO:
#            glDeleteBuffers(2, self.buffers)#

    cdef build_arrays(self):
        cdef double scale
        cdef vector3 V, N
        cdef GLfloat* NV
        cdef GLushort* T
#        glGenBuffers(2, self.buffers)
        NVsize = 6*len(self.mesh.vertices)*sizeof(GLfloat)
        self.nv_array = NV = <GLfloat *> malloc(NVsize)
        for vertex in self.mesh.vertices:
            scale = 1 + sqrt(max(0, 1 - vertex.norm_squared))
            V = vertex/scale
            N = self.center - V
            N = N/N.norm
            NV[0], NV[1], NV[2] = N.x, N.y, N.z
            NV[3], NV[4], NV[5] = V.x, V.y, V.z
            NV += 6
#        if self.useVBO:
#            glBindBuffer(GL_ARRAY_BUFFER, self.buffers[0])
#            glBufferData(GL_ARRAY_BUFFER, NVsize, self.nv_array,
#                         GL_STATIC_DRAW)
#            glBindBuffer(GL_ARRAY_BUFFER, 0)

        self.count = 3*len(self.mesh.triangles)
        Tsize = self.count*sizeof(GLushort)
        self.indices = T = <GLushort *> malloc(Tsize)
        for triangle in self.mesh.triangles:
            T[0], T[1], T[2] = triangle
            T += 3
#        if self.useVBO:
#            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self.buffers[1])
#            glBufferData(GL_ELEMENT_ARRAY_BUFFER, Tsize, self.indices,
#                         GL_STATIC_DRAW)
#            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0)

    def draw(self, use_material=True):
        if use_material:
            self.set_material()
#        if self.useVBO:
#            glBindBuffer(GL_ARRAY_BUFFER, self.buffers[0])
#            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self.buffers[1])
#            glNormalPointer(GL_FLOAT, 6*sizeof(GLfloat), <GLfloat*> NULL)
#            glVertexPointer(3, GL_FLOAT, 6*sizeof(GLfloat), <GLfloat*> NULL+3)
#        else:
#            glNormalPointer(GL_FLOAT, 6*sizeof(GLfloat), self.nv_array)
#            glVertexPointer(3, GL_FLOAT, 6*sizeof(GLfloat), self.nv_array+3)
        glNormalPointer(GL_FLOAT, 6*sizeof(GLfloat), self.nv_array)
        glVertexPointer(3, GL_FLOAT, 6*sizeof(GLfloat), self.nv_array+3)
        glEnableClientState(GL_NORMAL_ARRAY)
        glEnableClientState(GL_VERTEX_ARRAY)
#        if self.useVBO:
#            glDrawElements(GL_TRIANGLES, self.count, GL_UNSIGNED_SHORT,
#                           NULL)
#        else:
#            glDrawElements(GL_TRIANGLES, self.count, GL_UNSIGNED_SHORT,
#                           self.indices)
        glDrawElements(GL_TRIANGLES, self.count, GL_UNSIGNED_SHORT,
                       self.indices)
        glDisableClientState(GL_NORMAL_ARRAY)
        glDisableClientState(GL_VERTEX_ARRAY)
#        if self.useVBO:
#            glBindBuffer(GL_ARRAY_BUFFER, 0)
#            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0)

    def build_display_list(self, list_id):
        glNewList(list_id, GL_COMPILE) 
        self.draw()
        glEndList()

cdef class PoincarePolygon(GLobject):
    """
    Draws a geodesic polygon in the Poincare model. The geometric
    parameters are the vertex coordinates in the Klein model plus the
    coordinates of the center of the sphere which represents the plane
    of the triangle in the Poincare model.  The polygon is drawn by
    subdividing into Poincare Triangles by coning from the barycenter,
    then drawing each triangle.
    """
    cdef vertices, center, triangles

    def __init__(self, vertices, center, **kwargs):
        self.vertices = vertices
        self.center = center
        self.triangulate()

    def triangulate(self):
        Vlist = self.vertices
        zero = vector3((0,0,0))
        N = len(Vlist)
        self.triangles = []
        centroid = sum(Vlist, zero)/N
        for i in range(0,N):
            vertices = [centroid, Vlist[i-1],Vlist[i]]
            self.triangles.append(PoincareTriangle(vertices, self.center))

    def draw(self):
        self.set_material()
        for triangle in self.triangles:
            triangle.draw(use_material=False)

    def build_display_list(self, list_id):
        glNewList(list_id, GL_COMPILE) 
        self.draw()
        glEndList()

cdef class KleinPolygon(GLobject):
    """
    Draws a geodesic polygon in the Klein model. The geometric
    parameters are the vertex coordinates in the Klein model plus the
    coordinates of the nearest point to origin which lies on the plane
    containing the polygon.  The polygon is drawn as an OpenGL
    Polygon.
    """
    cdef vertices, closest

    def __init__(self, vertices, closest, **kwargs):
        self.vertices = vertices
        self.closest = closest

    def draw(self):
        N = self.closest/self.closest.norm
        self.set_material()
        glBegin(GL_POLYGON)
        glNormal3f(N.x, N.y, N.z)
        for V in self.vertices:
            glVertex3f(V.x, V.y, V.z)
        glEnd()

    def build_display_list(self, list_id):
        glNewList(list_id, GL_COMPILE) 
        self.draw()
        glEndList()

class HyperbolicPolyhedron:
   """
   A hyperbolic polyhedron for display in OpenGL, either in the Klein
   model or the Poincare model.  Includes a representation of the
   sphere at infinity.  It is initialized with the SnapPea description
   of the faces of a Dirichlet domain, represented as a list of
   dictionaries.
   """

   def __init__(self, facedicts, model_var, sphere_var):
     self.model = model_var
     self.sphere = sphere_var
     self.face_specular = [0.5, 0.5, 0.5, 1.0]
     self.front_shininess = 50.0
     self.back_shininess = 50.0
     self.sphere_listid = glGenLists(1)
     self.GLU = GLU_context()
     self.S_infinity = Sphere(GLU=self.GLU,
                              filled=False,
                              color=[1.0, 1.0, 1.0, .2],
                              front_specular=[0.5, 0.5, 0.5, 1.0],
                              front_shininess=50.0)
     self.S_infinity.build_display_list(self.sphere_listid, 1.0, 30, 30)
     self.Klein_faces = []
     self.Poincare_faces = []
     for dict in facedicts:
         vertices = [vector3(vertex) for vertex in dict['vertices']]
         closest = vector3(dict['closest'])
         center = closest*(1/dict['distance']**2)
         color = hls_to_rgb(dict['hue'], 0.5, 1.0) + (1.0,)
         self.Klein_faces.append(
             KleinPolygon(vertices, closest,
                             color=color,
                             front_specular=self.face_specular,
                             back_specular=self.face_specular,
                             front_shininess=self.front_shininess,
                             back_shininess=self.back_shininess))
         self.Poincare_faces.append(
             PoincarePolygon(vertices, center,
                             color=color,
                             front_specular=self.face_specular,
                             back_specular=self.face_specular,
                             front_shininess=self.front_shininess,
                             back_shininess=self.back_shininess))
     self.klein_list = glGenLists(1)
     self.build_klein_poly(self.klein_list)
     self.poincare_list = glGenLists(1)
     self.build_poincare_poly(self.poincare_list)

   def draw(self, *args):
       model = self.model.get()
       if model == 'Klein':
           glCallList(self.klein_list)
       elif model == 'Poincare':
           glCallList(self.poincare_list)
       if self.sphere.get():
           glCallList(self.sphere_listid)

   def build_klein_poly(self, list):
     glNewList(list, GL_COMPILE) 
     for face in self.Klein_faces:
       face.draw()
     glEndList()

   def build_poincare_poly(self, list):
     glNewList(list, GL_COMPILE) 
     for face in self.Poincare_faces:
       face.draw()
     glEndList()

class Colorizer:
    """
    Callable class which returns a color when passed an index.
    For indices larger than 5, a random color is generated and remembered.
    """
    basic_color = {0 : [1.0, 0.2, 0.2],
                   1 : [0.2, 0.2, 1.0],
                   2 : [0.2, 1.0, 0.2],
                   3 : [0.6, 0.6, 0.2],
                   4 : [0.2, 0.6, 0.6],
                   5 : [0.6, 0.6, 0.2]}

    def __init__(self, intensity=1.4, alpha=0.8):
        self.intensity = intensity
        self.alpha = alpha
        self.colors = {}

    def __call__(self, index):
        if 0 <= index < 6:
            R, G, B = self.basic_color[index]
        else:
            try:
                R, G, B = self.colors[index]
            except KeyError:
                R, G, B = random(), random(), random()
                scale = self.intensity/(R+G+B)
                self.colors[index] = [scale*R, scale*G, scale*B]
                R, G, B = self.colors[index]
        return [R, G, B, self.alpha]

GetColor = Colorizer()

cdef class HoroballGroup:
    """
    A fundamental set of horoballs for a single cusp.  The geometric
    parameter for the draw method is a list of shifts (M,L), meaning
    that each sphere should be drawn translated by M meridians and L
    longitudes.
    """
    cdef horoballs, spheres, meridian, longitude
    cdef GLUquadric* glu_quadric
    cdef GLfloat color[4]
    cdef double cutoff
    cdef int cusp_index
    
    def __init__(self, GLU_context GLU, horoballs, meridian, longitude):
        self.horoballs = horoballs
        self.meridian = meridian
        self.longitude = longitude
        self.glu_quadric = GLU.glu_quadric
        self.build_spheres()

    def build_spheres(self):
        self.spheres = []
        for D in self.horoballs:
            radius = D['radius']
            center = vector3((D['center'].real, D['center'].imag, D['radius']))
            color = GetColor(D['index'])
            self.spheres.append( (radius, center, color) )
            # Sort spheres by radius so smaller ones are drawn first.
            self.spheres.sort()

    def draw(self, shifts, full_list=True):
        for radius, center, color in self.spheres:
            for i from 0 <= i < 4:
                self.color[i] = color[i]
            glColor4fv(self.color)
            slices = max(20, int(60*radius))
            if full_list:
                stacks = slices
            else:
                stacks = 2
            glPushMatrix()
            glTranslatef(center.x, center.y, center.z)
            for M, L in shifts:
                disp = M*self.meridian + L*self.longitude
                glPushMatrix()
                glTranslatef(disp.real, disp.imag, 0.0)
                gluSphere(self.glu_quadric, radius, slices, stacks)
                glPopMatrix()
            glPopMatrix()

    def build_display_list(self, list_id, shifts, full_list=True):
        glNewList(list_id, GL_COMPILE) 
        self.draw(shifts, full_list)
        glEndList()

cdef class Parallelogram(GLobject):
    """
    Draws a parallelogram on the xy-plane centered at (0,0). The geometric
    parameters are complex numbers corresponding to the two side vectors.
    """

    def draw(self, s1, s2):
        glLineWidth(2.0)
        glColor4f(0.0, 0.0, 0.0, 1.0)
        glBegin(GL_LINE_LOOP)
        p = -(s1+s2)/2
        glVertex3f(p.real, p.imag, 0.0)
        p += s1
        glVertex3f(p.real, p.imag, 0.0)
        p += s2
        glVertex3f(p.real, p.imag, 0.0)
        p -= s1
        glVertex3f(p.real, p.imag, 0.0)
        glEnd()

    def build_display_list(self, list_id, s1, s2):
        glNewList(list_id, GL_COMPILE) 
        self.draw(s1, s2)
        glEndList()

cdef class FordEdgeSet:
    """
    A fundamental set of edges for the component of the Ford domain
    associated to a given cusp, projected to the xy-plane in upper
    half-space.  The geometric parameter for the draw method is a list
    of shifts (M,L), meaning that each segment should be drawn
    translated by M meridians and L longitudes.
    """
    cdef segments, longitude, meridian
    
    def __init__(self, segments, longitude, meridian):
        self.segments = segments 
        self.longitude, self.meridian = longitude, meridian

    def draw(self, shifts):
        glLineWidth(2.0)
        glColor4f(0.0, 0.0, 0.0, 1.0)
        for M, L in shifts:
            disp = M*self.meridian + L*self.longitude
            glPushMatrix()
            glTranslatef(disp.real, disp.imag, 0.0)
            for P1, P2 in self.segments:
                glBegin(GL_LINES)
                glVertex3f(P1.real, P1.imag, 0.0)
                glVertex3f(P2.real, P2.imag, 0.0)
                glEnd()
            glPopMatrix()

    def build_display_list(self, list_id, shifts):
        glNewList(list_id, GL_COMPILE) 
        self.draw(shifts)
        glEndList()

cdef class TriangulationEdgeSet:
    """
    A fundamental set of edges for the 1-skeleton of the canonical
    triangulation dual to the Ford domain, projected to the
    xy-plane in upper half-space.  The geometric parameter for the
    draw method is a list of shifts (M,L), meaning that each segment
    should be drawn translated by M meridians and L longitudes.
    """
    cdef segments, longitude, meridian
    
    def __init__(self, triangulation, longitude, meridian):
        self.segments = [D['endpoints'] for D in triangulation] 
        self.longitude, self.meridian = longitude, meridian

    def draw(self, shifts):
        glLineWidth(2.0)
        glColor4f(1.0, 1.0, 1.0, 1.0)
        for M, L in shifts:
            disp = M*self.meridian + L*self.longitude
            glPushMatrix()
            glTranslatef(disp.real, disp.imag, 0.0)
            for P1, P2 in self.segments:
                glBegin(GL_LINES)
                glVertex3f(P1.real, P1.imag, 0.0)
                glVertex3f(P2.real, P2.imag, 0.0)
                glEnd()
            glPopMatrix()

    def build_display_list(self, list_id, shifts):
        glNewList(list_id, GL_COMPILE) 
        self.draw(shifts)
        glEndList()

cdef class HoroballScene:
    """
    A family of Horoball Groups, one per cusp.  The variable which_cusp
    selects which group is visible.
    """
    cdef nbhd
    cdef meridian, longitude, offset
    cdef GLU, cusp_view, Ford, tri, pgram, shifts
    cdef pgram_var, Ford_var, tri_var
    cdef GLfloat Xangle, Yangle
    cdef GLint ball_list_id, pgram_list_id, Ford_list_id, tri_list_id
    cdef double cutoff
    cdef int which_cusp

    def __init__(self, nbhd, pgram_var, Ford_var, tri_var,
                 cutoff=0.1, which_cusp=0):
        self.nbhd = nbhd
        self.tri_var = tri_var
        self.Ford_var = Ford_var
        self.which_cusp = which_cusp
        self.pgram_var = pgram_var
        self.offset = 0.0j
        self.Xangle, self.Yangle = 0.0, 0.0
        self.GLU = GLU_context()
        self.setup_quadric(self.GLU)
        self.pgram = Parallelogram()
        self.set_cutoff(cutoff)
        self.build_scene()

    def set_cutoff(self, cutoff):
        self.cutoff = cutoff
        
    def build_scene(self, full_list=True):
        self.meridian, self.longitude = self.nbhd.translations(
            self.which_cusp)
        self.build_shifts()
        self.cusp_view = HoroballGroup(
            self.GLU,
            self.nbhd.horoballs(
                self.cutoff,
                self.which_cusp,
                full_list),
            self.meridian,
            self.longitude)
        self.Ford = FordEdgeSet(
                self.nbhd.Ford_domain(self.which_cusp),
                self.longitude, self.meridian)
        self.tri = TriangulationEdgeSet(
                self.nbhd.triangulation(self.which_cusp),
                self.longitude, self.meridian)
        self.gl_compile(full_list)

    def setup_quadric(self, GLU_context GLU):
        gluQuadricDrawStyle(GLU.glu_quadric, GLU_FILL)
        gluQuadricNormals(GLU.glu_quadric, GLU_SMOOTH)

    def build_shifts(self):
        size = 2.1*max(self.longitude.real, self.meridian.imag)
        M = int(ceil(size/abs(self.meridian.imag)))
        N = int(ceil(size/self.longitude.real))
        self.shifts = []
        for m in range(-M,M):
            shear = m*self.meridian.real/self.longitude.real
            left = int(floor(-shear-N))
            for n in range(left,left+2*N):
                self.shifts.append((m,n))
    
    cdef change_basis(self, z):
        cdef GLdouble model[16]
        cdef GLdouble A, B, C, D, Z
        glGetDoublev(GL_MODELVIEW_MATRIX, model)
        Z = 1/(model[12] + model[15])
        A = (model[0] + model[3])*Z
        C = (model[4] + model[7])*Z
        Z = 1/(model[13] + model[15])
        B = (model[1] + model[3])*Z
        D = (model[5] + model[7])*Z
        w = (D*z.real - C*z.imag) + (-B*z.real + A*z.imag)*1j 
        return w*abs(z)/abs(w)

    def translate(self, z):
        """
        Translate modulo the cusp stabilizer.
        """
        z = self.change_basis(z)
        z += self.offset
        z = z - (z.imag//self.meridian.imag)*self.meridian
        z = z - (z.real//self.longitude.real)*self.longitude
        self.offset = z

    def get_lists(self):
        if self.pgram_list_id != 0:
            glDeleteLists(self.pgram_list_id, 4)
        self.pgram_list_id = glGenLists(4)
        self.ball_list_id = self.pgram_list_id + 1
        self.Ford_list_id = self.pgram_list_id + 2
        self.tri_list_id = self.pgram_list_id + 3

    def gl_compile(self, full_list):
        self.get_lists()
        self.pgram.build_display_list(self.pgram_list_id,
                                      self.longitude, self.meridian)
        self.cusp_view.build_display_list(self.ball_list_id,
                                          self.shifts, full_list)
        self.Ford.build_display_list(self.Ford_list_id, self.shifts)
        self.tri.build_display_list(self.tri_list_id, self.shifts)

    def draw_segments(self):
        glPushMatrix()
        glTranslatef(self.offset.real, self.offset.imag, 0.0)
        if self.Ford_var.get():
            glCallList(self.Ford_list_id)
        if self.tri_var.get():
            glCallList(self.tri_list_id)
        glPopMatrix()
        if self.pgram_var.get():
            glCallList(self.pgram_list_id)


    def draw(self, *args):
        """
        The scene is drawn translated by self.offset, but the
        parallelogram stays fixed.
        """
        glPushMatrix()
        glTranslatef(self.offset.real, self.offset.imag, 0.0)
        glCallList(self.ball_list_id)
        glPopMatrix()
        # Draw segments without a depth buffer, for better anti-aliasing
        glDisable(GL_DEPTH_TEST)
        self.draw_segments()
        glEnable(GL_DEPTH_TEST)
        # Do it again to make them a little darker
        self.draw_segments()

# Methods to translate and rotate our scene.

cdef glTranslateScene(s, x, y, mousex, mousey):
    cdef GLfloat X, Y
    cdef GLdouble mat[16]

    X, Y = s * (x - mousex), s * (mousey - y)
#    glMatrixMode(GL_MODELVIEW)
    glGetDoublev(GL_MODELVIEW_MATRIX, mat)
    glLoadIdentity()
    glTranslatef(X, Y, 0.0)
    glMultMatrixd(mat)

cdef glRotateScene(xcenter, ycenter, zcenter, Xangle, Yangle, cos_bound=-1.0):
    cdef GLdouble mat[16]

    #glMatrixMode(GL_MODELVIEW)
    glGetDoublev(GL_MODELVIEW_MATRIX, mat)
    glLoadIdentity()
    glTranslatef(xcenter, ycenter, zcenter)
    glRotatef(Yangle, 1., 0., 0.)
    glRotatef(Xangle, 0., 1., 0.)
    glTranslatef(-xcenter, -ycenter, -zcenter)
    glMultMatrixd(mat)
    if check_angle(cos_bound) is False:
        glLoadMatrixd(mat)

cdef check_angle(cos_bound):
    cdef GLdouble mat[16]

    glGetDoublev(GL_MODELVIEW_MATRIX, mat)
    cosine = (mat[10] + mat[14])/(mat[11] + mat[15])
    if abs(cosine) > cos_bound:
        return True
    return False

class RawOpenGLWidget(Tk_.Widget, Tk_.Misc):
    """
    Widget without any sophisticated bindings
    by Tom Schwaller
    """

    def __init__(self, master, cnf={}, **kw):
        snappy_dir = os.path.dirname(__file__)
        # Hack to make py2exe behave:
        if not snappy_dir.endswith('snappy'):
            snappy_dir = os.path.join(snappy_dir, 'snappy')

        curr_platform = sys.platform
        if curr_platform[:5] == "linux" and platform.architecture()[0] == '64bit':
            curr_platform += "-x86_64"
        Togl_path = os.path.join( snappy_dir,
                              curr_platform + "-tk" + master.getvar("tk_version"))
        master.tk.call('lappend', 'auto_path', Togl_path)
        master.tk.call('package', 'require', 'Togl')

        Tk_.Widget.__init__(self, master, 'togl', cnf, kw)
        self.root = master
        self.bind('<Map>', self.tkMap)
        self.bind('<Expose>', self.tkExpose)
        self.bind('<Configure>', self.tkExpose)

    def tkRedraw(self, *dummy):
        self.tk.call(self._w, 'makecurrent')
        glPushMatrix()
        self.update_idletasks()
        self.redraw()
        glFlush()
        glPopMatrix()
        self.tk.call(self._w, 'swapbuffers')

    def tkMap(self, *dummy):
        self.tkExpose()

    def tkExpose(self, *dummy):
        self.tkRedraw()

class OpenGLWidget(RawOpenGLWidget):
    """
    Tkinter bindings for an OpenGL widget.
    Mike Hartshorn
    Department of Chemistry
    University of York, UK
    http://www.yorvic.york.ac.uk/~mjh/
    """

    def __init__(self, master=None, help='No help is available.',
                 mouse_pick=False, mouse_rotate=True, mouse_translate=False,
                 mouse_scale=False, translate=None, cos_bound=-1.0,
                 cnf={}, **kw):
        """
        Create an opengl widget.  Arrange for redraws when the window is
        exposed or when it changes size.
        """

        apply(RawOpenGLWidget.__init__, (self, master, cnf), kw)
        if translate:
            self.translate = translate
        else:
            self.translate = self.tkTranslate
        self.cos_bound = cos_bound
        self.help_text = help
        self.initialised = 0
        if sys.platform == 'darwin':
            self.config(cursor='hand')
        else:
            self.config(cursor='fleur')

        # Current coordinates of the mouse.
        self.xmouse = 0
        self.ymouse = 0

        # Where we are centering.
        self.xcenter = 0.0
        self.ycenter = 0.0
        self.zcenter = 0.0

        # The _back color
        self.r_back = 1.
        self.g_back = 0.
        self.b_back = 1.

        # Where the eye is
        self.distance = 10.0
    
        # Field of view in y direction
        self.fovy = 30.0

        # Position of clipping planes.
        self.near = 1.0
        self.far = 100.0

        # Is the widget allowed to autospin?
        self.autospin_allowed = 0

        # Is the widget currently autospinning?
        self.autospin = 0

        # Dictionary of key actions (keysym:function) .
        self.key_action = {}
        
        # Bindings for events.
        self.bind('<Map>', self.tkMap)
        self.bind('<Expose>', self.tkExpose)
        self.bind('<Configure>', self.tkExpose)
        if mouse_pick:
            self.bind('<Control-Button-1>', self.tkHandlePick)
            self.bind('<Control-Button-1><ButtonRelease-1>', self.tkHandlePick)
        if mouse_translate and mouse_rotate:
            self.bind('<Button-1>', self.tkRecordMouse)
            self.bind('<B1-Motion>', self.translate)
            if sys.platform == 'darwin':
                self.bind('<Shift-Button-1>', self.StartRotate)
                self.bind('<Shift-B1-Motion>', self.tkRotate)
                self.bind('<ButtonRelease-1>', self.tkAutoSpin)
            else:
                self.bind('<Button-3>', self.StartRotate)
                self.bind('<B3-Motion>', self.tkRotate)
                self.bind('<ButtonRelease-3>', self.tkAutoSpin)
        elif mouse_rotate:
            self.bind('<Button-1>', self.StartRotate)
            self.bind('<B1-Motion>', self.tkRotate)
            self.bind('<ButtonRelease-1>', self.tkAutoSpin)
        elif mouse_translate:
            self.bind('<Button-1>', self.tkRecordMouse)
            self.bind('<B1-Motion>', self.translate)
        if mouse_scale:
            self.bind('<Button-2>', self.tkRecordMouse)
            self.bind('<B2-Motion>', self.tkScale)
            self.bind('<KeyPress>', self.tkKeyPress)

    def help(self):
        """
        Help message for the widget.
        """
        tkMessageBox.showinfo('Viewer Help', self.help_text)

    def activate(self):
        """
        Cause this OpenGLWidget to be the current destination for
        drawing.  Does NOT make the widget be the focus of keyboard
        events; SnapPy OpenGL widgets to not accept keyboard events.
        """
        self.tk.call(self._w, 'makecurrent')
        #self.focus_set()

    def set_background(self, r, g, b):
        """
        Change the background colour of the widget.
        """
        self.r_back = r
        self.g_back = g
        self.b_back = b
        self.tkRedraw()

    def set_centerpoint(self, x, y, z):
        """
        Set the new center point for the model.
        This is where we are looking.
        """
        self.xcenter = x
        self.ycenter = y
        self.zcenter = z
        self.tkRedraw()

    def set_eyepoint(self, distance):
        """
        Set how far the eye is from the position we are looking.
        """
        self.distance = distance
        self.tkRedraw()

    def reset(self, redraw=True):
        """
        Reset rotation matrix for this widget.
        """
        self.autospin = 0
        self.activate()
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity()
        if redraw:
            self.tkRedraw()

    def tkHandlePick(self, event):
        """
        Handle a pick on the scene.
        """
        cdef GLdouble objX, objY, objZ
        cdef GLdouble model[16], proj[16]
        cdef GLint view[4]

        if hasattr(self, 'pick'):
          # here we need to use glu.UnProject
          # Tk and X have their origin top left, 
          # while OpenGLWidget has its origin bottom left.
          # So we need to subtract y from the window height to get
          # the proper pick position for OpenGLWidget
            realy = self.winfo_height() - event.y
            glGetDoublev(GL_MODELVIEW_MATRIX, model)
            glGetDoublev(GL_PROJECTION_MATRIX, proj)
            glGetIntegerv(GL_VIEWPORT, view)
            gluUnProject(event.x, realy, 0., model, proj, view, &objX, &objY, &objZ)
            p1 = (objX, objY, objZ)
            gluUnProject(event.x, realy, 1., model, proj, view, &objX, &objY, &objZ)
            p2 = (objX, objY, objZ)

            if self.pick(self, p1, p2):
                # If the pick method returns true we redraw the scene.
                self.tkRedraw()

    def tkRecordMouse(self, event):
        """
        Record the current mouse position.
        """
        self.xmouse = event.x
        self.ymouse = event.y

    def StartRotate(self, event):
        # Switch off any autospinning if it was happening
        self.autospin = 0
        self.tkRecordMouse(event)

    def tkScale(self, event):
        """
        Scale the scene.  Achieved by moving the eye position.
        """
        scale = 1 - 0.01 * (event.y - self.ymouse)
        self.distance = self.distance * scale
        self.tkRedraw()
        self.tkRecordMouse(event)

    def do_AutoSpin(self):
        s = 0.1
        self.activate()

        glRotateScene(self.xcenter, self.ycenter, self.zcenter,
                      s*self.yspin, s*self.xspin)
        self.tkRedraw()

        if self.autospin:
            self.after(10, self.do_AutoSpin)

    def tkAutoSpin(self, event):
        """
        Perform autospin of scene.
        """
        self.after(16)
        self.update_idletasks()
        x = self.tk.getint(self.tk.call('winfo', 'pointerx', self._w))
        y = self.tk.getint(self.tk.call('winfo', 'pointery', self._w))

        if self.autospin_allowed:
            if x != event.x_root and y != event.y_root:
                self.autospin = 1

            self.yspin = x - event.x_root
            self.xspin = y - event.y_root
            self.after(10, self.do_AutoSpin)

    def tkRotate(self, event):
        """
        Perform rotation of scene.
        """
        cdef GLfloat Xangle, Yangle
        self.activate()
        Xangle = 0.5 * (event.x - self.xmouse)
        Yangle = 0.5 * (event.y - self.ymouse)
        glRotateScene(self.xcenter, self.ycenter, self.zcenter,
                      Xangle, Yangle, self.cos_bound)
        self.tkRedraw()
        self.tkRecordMouse(event)

    def tkTranslate(self, event):
        """
        Perform translation of scene.
        """
        self.activate()
        glTranslateScene(0.05, event.x, event.y, self.xmouse, self.ymouse)
        self.tkRedraw()
        self.tkRecordMouse(event)

    def flip(self):
        """
        Rotate 180 degrees about the x-axix.
        """
        self.activate()
        glRotateScene(self.xcenter, self.ycenter, self.zcenter,
                      180.0, 0.0)
        self.tkRedraw()

    def mouse_update(self, event):
        """
        Redraw the scene and save the mouse coordinates.
        """
        self.tkRedraw()
        self.tkRecordMouse(event)

    def tkRedraw(self, *dummy):
        """Cause the opengl widget to redraw itself."""
        if not self.initialised: return
        self.activate()
        glPushMatrix()                        # Protect our matrix
        self.update_idletasks()
        w = self.winfo_width()
        h = self.winfo_height()
        glViewport(0, 0, w, h)

        # Clear the background and depth buffer.
        glClearColor(self.r_back, self.g_back, self.b_back, 0.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity()
        gluPerspective(self.fovy, float(w)/float(h), self.near, self.far)
        gluLookAt(self.xcenter, self.ycenter, self.zcenter + self.distance,
                  self.xcenter, self.ycenter, self.zcenter, 0.0, 1.0, 0.0)
        glMatrixMode(GL_MODELVIEW);

        # Call objects redraw method.
        try:
            self.redraw(self)
        except AttributeError:
            pass
        glFlush()                                # Tidy up
        glPopMatrix()                            # Restore the matrix

        self.tk.call(self._w, 'swapbuffers')

    def tkMap(self, *dummy):
        """
        Cause the opengl widget to redraw itself.
        """
        self.tkExpose()

    def tkExpose(self, *dummy):
        """
        Redraw the widget.  Make it active, update tk events, call redraw
        procedure and swap the buffers.  Note: swapbuffers is clever
        enough to only swap double buffered visuals.
        """
        self.activate()
        if not self.initialised:
            self.initialised = 1
        self.tkRedraw()

    def tkKeyPress(self, event):
        """
        Handle keyboard events.
        """
        try:
            self.key_action[event.keysym]()
        except KeyError:
            pass
        if not self.autospin:
            self.tkRedraw()

    def tkPrint(self, file):
        """
        Turn the current scene into PostScript via the feedback buffer.
        """
        self.activate()
