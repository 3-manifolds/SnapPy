include "CyOpenGL.pxi"

# This is part of the UCS2 hack.
cdef public UCS2_hack (char *string, Py_ssize_t length, char *errors) :   
    return string 

from .infodialog import InfoDialog
from . import togl
    
import os, sys, platform
from colorsys import hls_to_rgb
from math import sqrt, ceil, floor, pi, sin, cos, tan
from random import random

if sys.version_info[0] < 3: 
    import Tkinter as Tk_
else:
    import tkinter as Tk_

def GetString(string):
    enumdict = {
        'GL_VENDOR': GL_VENDOR,
        'GL_RENDERER': GL_RENDERER,
        'GL_VERSION': GL_VERSION,
        'GL_EXTENSIONS': GL_EXTENSIONS}
    try:
        result = <const char*>glGetString(enumdict[string])
    except KeyError:
        raise ValueError(
            "Invalid enum. Must be 'GL_VENDOR', 'GL_RENDERER', 'GL_VERSION', or 'GL_EXTENSIONS'")
    if result:
        return result
    else:
        raise RuntimeError('No result - is there a current OpenGL context?')

cdef class vector3:
    """
    A simple real 3-dimensional vector which supports addition,
    subtraction and right multiplication or division by scalars.
    Attributes include its norm and the square of its norm.
    """
    cdef readonly double x, y, z, norm_squared, norm

    def __cinit__(self, triple):
        self.x, self.y, self.z = triple
        self.norm_squared = self.x*self.x + self.y*self.y + self.z*self.z 
        self.norm = sqrt(self.norm_squared)

    def __repr__(self):
        """
        >>> vector3( (0, 1, 2) )
        < 0.0, 1.0, 2.0 >
        """
        return '< %s, %s, %s >'%(self.x, self.y, self.z)

    def __add__(self, vector):
        return vector3([self.x+vector.x, self.y+vector.y, self.z+vector.z])

    def __sub__(self, vector):
        return vector3([self.x-vector.x, self.y-vector.y, self.z-vector.z])

    def __mul__(self, scalar):
        return vector3([self.x*scalar, self.y*scalar, self.z*scalar])

    def __div__(self, scalar):
        return vector3([self.x/scalar, self.y/scalar, self.z/scalar])

    def __truediv__(self, scalar):
        return vector3([self.x/scalar, self.y/scalar, self.z/scalar])

cdef class GL_context:
    """
    Sets up our default OpenGL environment.
    """

    def __cinit__(self):
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
    cdef GLfloat emission[4]
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
            self.color[n] = color[n]
            self.front_specular[n] = front_specular[n]
            self.back_specular[n] = back_specular[n]
            self.emission[n] = 0.0
        self.front_shininess = front_shininess
        self.back_shininess = back_shininess

    def set_material(self):
        glMaterialfv(GL_FRONT, GL_SPECULAR, self.front_specular)
        glMaterialf(GL_FRONT, GL_SHININESS, self.front_shininess)
        glMaterialfv(GL_BACK,  GL_SPECULAR, self.back_specular)
        glMaterialf(GL_BACK,  GL_SHININESS, self.back_shininess)
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, self.emission)
        glColor4fv(self.color)

    def draw(self, *args, **kwargs):
        """
        Issue the OpenGL commands to draw this object.
        (Subclasses must override this.)
        """

    def build_display_list(self, list_id, *args, **kwargs):
        """
        Generate a display list containing the commands to draw this object.
        """
        glNewList(list_id, GL_COMPILE) 
        self.draw(*args, **kwargs)
        glEndList()

cdef class WireframeSphere(GLobject):
    """
    Draw a wireframe sphere.
    """

    def draw(self, GLfloat radius, GLint slices, GLint stacks):
        assert slices % 2 == 0 and stacks % 2 == 0
        self.set_material()
        r = radius
        # We put the north pole on the y-axis. 
        glPushMatrix()
        glLoadIdentity()
        glRotatef(90, 1.0, 0.0, 0.0)

        dtheta = 2*pi/slices
        dphi = pi/stacks
        def draw_point(phi, theta):
            x, y, z = sin(phi)*sin(theta), sin(phi)*cos(theta), cos(phi)
            glNormal3f(x, y, z)
            glVertex3f(r*x, r*y, r*z)
                
        # Draw the longitudes
        for i in range(slices):
            theta = dtheta*i
            glBegin(GL_LINE_STRIP)
            for j in range(stacks + 1):
                draw_point(dphi*j, theta)
            glEnd()

        # Draw the latitudes
        for j in range(1, stacks):
            glBegin(GL_LINE_LOOP)
            phi = dphi*j
            for i in range(0, slices):
                draw_point(phi, dtheta*i)
            glEnd()
        glPopMatrix()   

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


cdef class MeshedSurface(GLobject):
    """ 
    An object made out of a triangular mesh. See the subclass
    Horosphere below for a typical example.
    """
    
    cdef vertices, normals, triangles, count
    cdef GLfloat* nv_array
    cdef GLushort* indices

    def __dealloc__(self):
        free(self.nv_array)
        free(self.indices)

    cdef build_arrays(self):
        cdef double scale
        cdef vector3 V, N
        cdef GLfloat* NV
        cdef GLushort* T
        NVsize = 6*len(self.vertices)*sizeof(GLfloat)
        self.nv_array = NV = <GLfloat *> malloc(NVsize)
        for V, N in zip(self.vertices, self.normals):
            NV[0], NV[1], NV[2] = N.x, N.y, N.z
            NV[3], NV[4], NV[5] = V.x, V.y, V.z
            NV += 6

        self.count = 3*len(self.triangles)
        Tsize = self.count*sizeof(GLushort)
        self.indices = T = <GLushort *> malloc(Tsize)
        for triangle in self.triangles:
            T[0], T[1], T[2] = triangle
            T += 3
    
    def draw(self, use_material=True):
        glNormalPointer(GL_FLOAT, 6*sizeof(GLfloat), self.nv_array)
        glVertexPointer(3, GL_FLOAT, 6*sizeof(GLfloat), self.nv_array+3)
        glEnableClientState(GL_NORMAL_ARRAY)
        glEnableClientState(GL_VERTEX_ARRAY)
        glDisableClientState(GL_COLOR_ARRAY)
        if use_material:
            self.set_material()
        glDrawElements(GL_TRIANGLES, self.count, GL_UNSIGNED_SHORT,
                       self.indices)
        glDisableClientState(GL_NORMAL_ARRAY)
        glDisableClientState(GL_VERTEX_ARRAY)
    
    

cdef class PoincareTriangle(MeshedSurface):
    """
    Draws a geodesic triangle in the Poincare model.  The geometric
    parameters are the vertex coordinates in the Klein model plus the
    coordinates of the center of the sphere which represents the plane
    of the triangle in the Poincare model.  The Poincare vertices are
    constructed by projecting the Klein vertices onto the sphere from
    the center.

    The triangle is drawn as a mesh, using vertex and index arrays.
    """
    cdef center, mesh, original_vertices

    def __init__(self, vertices, center, subdivision_depth=4, **kwargs):
        self.original_vertices = vertices
        self.center = center
        self.mesh = TriangleMesh(vertices)
        for n in range(subdivision_depth):
            self.mesh.subdivide()


        self.triangles = self.mesh.triangles
        self.vertices = vertices = []
        self.normals = normals = []
        for vertex in self.mesh.vertices:
            scale = 1 + sqrt(max(0, 1 - vertex.norm_squared))
            V = vertex/scale
            N = self.center - V
            N = N/N.norm
            vertices.append(V)
            normals.append(N)

        self.build_arrays()

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

class HyperbolicPolyhedron:
   """
   A hyperbolic polyhedron for display in OpenGL, either in the Klein
   model or the Poincare model.  Includes a representation of the
   sphere at infinity.  It is initialized with the SnapPea description
   of the faces of a Dirichlet domain, represented as a list of
   dictionaries.
   """

   def __init__(self, facedicts, model_var, sphere_var):
     self.facedicts = facedicts
     self.model = model_var
     self.sphere = sphere_var
     self.face_specular = [0.5, 0.5, 0.5, 1.0]
     self.front_shininess = 50.0
     self.back_shininess = 50.0
     self.sphere_list_id = glGenLists(1)
     self.S_infinity = WireframeSphere(color=[1.0, 1.0, 1.0, .2],
                              front_specular=[0.5, 0.5, 0.5, 1.0],
                              front_shininess=50.0)
     self.S_infinity.build_display_list(self.sphere_list_id, 1.0, 30, 30)
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
     self.klein_list_id = glGenLists(1)
     self.build_klein_poly(self.klein_list_id)
     self.poincare_list_id = glGenLists(1)
     self.build_poincare_poly(self.poincare_list_id)

   def destroy(self):
       glDeleteLists(self.sphere_list_id, 1)
       glDeleteLists(self.klein_list_id, 1)
       glDeleteLists(self.poincare_list_id, 1)
       
   def draw(self, *args):
       model = self.model.get()
       if model == 'Klein':
           glCallList(self.klein_list_id)
       elif model == 'Poincare':
           glCallList(self.poincare_list_id)
       if self.sphere.get():
           glCallList(self.sphere_list_id)

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

cdef class Colorizer:
    """
    Callable class which returns a color when passed an index.
    Uses the same algorithm as the SnapPea kernel.
    """
    cdef int base_hue[6]
    cdef double lightness, saturation, alpha

    def __cinit__(self):
        cdef int n
        # red blue green cyan magenta yellow
        cdef hues = [0,4,2,3,5,1]
        # maybe one day Cython will let you initialize C arrays
        for n in range(6):
            self.base_hue[n] = hues[n]

    def __init__(self, lightness=0.6, saturation=0.9, alpha=0.8):
        self.lightness = lightness
        self.saturation = saturation
        self.alpha = alpha

    def __call__(self, index):
        cdef double hue, R, G, B
        hue = (self.base_hue[index%6] + self.index_to_hue(index//6)) / 6.0
        R, G, B = self.hls_to_rgb(hue, self.lightness, self.saturation)  
        return [R, G, B, self.alpha]

    cdef double index_to_hue(self, int index):
        cdef unsigned int num=0, den=1
        while index:
            num = num<<1
            den = den<<1
            if index & 0x1:
                num += 1
            index = index>>1
        return <double>num/<double>den

    cdef hls_to_rgb(self, double h, double l, double s):
        if s == 0.0:
            return l, l, l
        if l <= 0.5:
            m2 = l * (1.0+s)
        else:
            m2 = l+s-(l*s)
        m1 = 2.0*l - m2
        return (self.hls_interp(m1, m2, h+1.0/3.0),
                self.hls_interp(m1, m2, h),
                self.hls_interp(m1, m2, h-1.0/3.0))

    cdef hls_interp(self, double m1, double m2, double hue):
        hue = hue % 1.0
        if hue < 1.0/6.0:
            return m1 + (m2-m1)*hue*6.0
        if hue < 0.5:
            return m2
        if hue < 2.0/3.0:
            return m1 + (m2-m1)*(2.0/3.0-hue)*6.0
        return m1

GetColor = Colorizer()

cdef class Horosphere(MeshedSurface):
    """
    Draw a horosphere.
    """
    cdef GLdouble radius
    cdef GLint stacks, slices

    def __init__(self,
                 color=[0.8,0.0,0.0,0.3],
                 radius=1.0,
                 front_specular = [0.8, 0.8, 0.8, 1.0], 
                 back_specular = [0.8, 0.8, 0.8, 1.0],
                 front_shininess = 50.0,
                 back_shininess = 0.0
                 ):
        self.radius = radius
        self.stacks = 2*max(2, int(8*radius))
        self.slices = max(20, int(60*radius))
        self.build_vertices_and_normals()
        self.build_triangles()
        self.build_arrays()

    def build_vertices_and_normals(self):
        a, b = self.stacks, self.slices
        dtheta = 2*pi/b
        dphi = pi/a

        verts = []
        phi = 0
        for j in range(a - 1):
            phi += dphi
            r, z = sin(phi), cos(phi)
            theta = 0
            for i in range(0, b):
                verts.append((r*cos(theta), r*sin(theta), z))
                theta += dtheta
        verts += [(0, 0, 1), (0, 0, -1)]
        assert len(verts) == a*b - b + 2

        self.vertices, self.normals = [], []
        for x, y, z in verts:
            N = vector3((x, y, z))
            V = N*self.radius
            self.vertices.append(V)
            self.normals.append(N)
            
    def build_triangles(self):
        self.triangles = tri = []
        a, b = self.stacks, self.slices

        # Start with the two polar caps
        north = a*b - b
        south = north + 1
        tri += [(north, i, (i + 1) % b) for i in range(b)]
        shift = north - b
        tri += [(south, shift + (i + 1) % b, shift + i) for i in range(b)]
        
        # Now build the rest with annular bands
        annulus = []
        for v0 in range(0, b):
            w0 = (v0 + 1) % b
            v1, w1 = v0 + b, w0 + b
            annulus += [(v0, v1, w1), (w0, v0, w1)]

        for s in range(0, b*(a - 2), b):
            tri += [(u + s, v + s, w + s) for u, v, w in annulus]
        
        
cdef class HoroballGroup:
    """
    A fundamental set of horoballs for a single cusp.  The paremeters
    R and T for the draw method are the coordinates of the right top
    corner of the visible rectangle with margins.  For each horoball,
    all meridian and longitude translations centered in the rectangle
    are drawn.
    """
    cdef horoballs, meridian, longitude,
    cdef keys, centers, spheres, list_ids
    cdef GLfloat color[4]
    cdef GLuint list_id_base, num_lists
    cdef double cutoff
    cdef original_indices
    
    def __init__(self, horoballs, indices, meridian, longitude):
        self.horoballs = horoballs
        self.meridian = complex(meridian)
        self.longitude = complex(longitude)
        self.list_id_base = 0
        self.num_lists = 0
        self.original_indices = indices
        self.build_spheres()

    def get_list_ids(self, N):
        self.delete_lists()
        self.list_id_base = glGenLists(N)
        self.num_lists = N

    def delete_lists(self):
        if self.list_id_base > 0:
            glDeleteLists(self.list_id_base, self.num_lists)
        
    def build_spheres(self):
        self.keys = keys = []
        self.spheres = spheres = {}
        self.centers = centers = {}
        self.list_ids = list_ids = {}
        for D in self.horoballs:
            z_center = D['center']
            radius = round(D['radius'], 10)
            index = D['index']
            key = (radius, index)
            center = vector3((z_center.real, z_center.imag, radius))
            color = GetColor(self.original_indices[index])
            try:
                centers[key].append(center)
            except KeyError:
                keys.append(key)
                centers[key] = [center]
                spheres[key] = Horosphere(radius=radius, color=color)
        keys.sort()
        self.get_list_ids(len(keys))
        n = self.list_id_base
        for key in keys:
            spheres[key].build_display_list(n)
            list_ids[key] = n
            n += 1

    def draw(self, R, T):
        vx, vy = self.meridian.real, self.meridian.imag
        ux = self.longitude.real
        for key in self.keys:
            list_id = self.list_ids[key]
            for center in self.centers[key]:
                x, y = center.x, center.y
                glPushMatrix()
                glTranslatef(x, y, center.z)
                N_min = -ceil( (T + y)/vy )
                N_max = ceil( (T - y)/vy )
                for n from N_min <= n <= N_max:
                    xn = x + n*vx
                    yn = y + n*vy
                    M_min = -ceil( (R + xn)/ux )
                    M_max = ceil( (R - xn)/ux )
                    for m from M_min <= m <= M_max:
                        disp = n*self.meridian + m*self.longitude
                        glPushMatrix()
                        glTranslatef(disp.real, disp.imag, 0.0)
                        glCallList(list_id)
                        glPopMatrix()
                glPopMatrix()

    def build_display_list(self, list_id, R, T):
        glNewList(list_id, GL_COMPILE) 
        self.draw(R, T)
        glEndList()

cdef class Parallelogram(GLobject):
    """
    Draws a parallelogram on the xy-plane centered at (0,0). The geometric
    parameters are complex numbers corresponding to the two side vectors.
    """

    def draw(self, s1, s2):
        glDisable(GL_LIGHTING)
        glLineWidth(2.0)
        glColor4f(1.0, 0.0, 1.0, 1.0)
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
        glEnable(GL_LIGHTING)

cdef class EdgeSet(GLobject):
    """
    Base class for collections of segments to be drawn into a
    horoball scene.

    The geometric parameter for the draw method is a list
    of shifts (M,L), meaning that each segment should be drawn
    translated by M meridians and L longitudes.
    """
    cdef segments, longitude, meridian, stipple
    
    def __init__(self, segments, longitude, meridian):
        self.segments = segments
        self.longitude, self.meridian = complex(longitude), complex(meridian)
        self.stipple = True

    cdef set_dark_color(self):
        glColor4f(0.0, 0.0, 0.0, 1.0)

    cdef set_light_color(self):
        glColor4f(0.7, 0.7, 0.7, 1.0)
    
    def draw(self, shifts, dark=True):
        glDisable(GL_LIGHTING)
        glLineWidth(2.0)
        if self.stipple:
            glEnable(GL_LINE_STIPPLE)
            glLineStipple(1, 0xcccc)
        if dark:
            self.set_dark_color()
        else:
            self.set_light_color()
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
        if self.stipple:
            glDisable(GL_LINE_STIPPLE)
        glEnable(GL_LIGHTING)

cdef class FordEdgeSet(EdgeSet):
    """
    A fundamental set of edges for the component of the Ford domain
    associated to a given cusp, projected to the xy-plane in upper
    half-space.
    """

cdef class TriangulationEdgeSet(EdgeSet):
    """
    A fundamental set of edges for the 1-skeleton of the canonical
    triangulation dual to the Ford domain, projected to the
    xy-plane in upper half-space.
    """
    def __init__(self, triangulation, longitude, meridian):
        self.segments = [D['endpoints'] for D in triangulation] 
        self.longitude, self.meridian = complex(longitude), complex(meridian)
        self.stipple = False

    cdef set_light_color(self):
        glColor4f(0.6, 0.6, 0.6, 1.0)

cdef class Label:
    cdef GLfloat x, y
    cdef codes
    
    def __init__(self, position, int_label):
        self.x, self.y = position.real, position.imag
        self.codes = [ord(c) for c in repr(int_label)]

    cdef get_shape(self):
        cdef SnapPy_glyph* glyph
        width, height = 0, 0
        for c in self.codes:
            glyph = SnapPy_font[c]
            width += glyph.width
            height = glyph.height
        return width, height

    def draw(self):
        cdef SnapPy_glyph* glyph
        glRasterPos2f(self.x, self.y)
        width, height = self.get_shape()
        # This is a trick to move the raster position in units of 1 pixel
        glBitmap(0, 0, 0, 0, -width/2, -height/2, NULL)
        for c in self.codes:
            glyph = SnapPy_font[c]
            if glyph != NULL:
                glDrawPixels(glyph.width, glyph.height,
                             GL_RGBA, GL_UNSIGNED_BYTE,
                             <GLvoid*> glyph.pixel_data)
                glBitmap(0, 0, 0, 0, glyph.width, 0, NULL)
    
cdef class LabelSet(GLobject):
    """
    Renders edge and vertex labels in the SnapPy font.
    """
    cdef segments, vertices, longitude, meridian, codes
    cdef SnapPy_glyph* glyph
    cdef GLfloat pix, x, y
    cdef int width, height

    def __init__(self, triangulation, longitude, meridian):
        self.longitude, self.meridian = complex(longitude), complex(meridian)
        self.segments = [ Label(sum(D['endpoints'])/2, D['indices'][1])
                          for D in triangulation]

        vertices = [ (complex(D['endpoints'][0]), int(D['indices'][0]))
                     for D in triangulation]
        vertices += [ (complex(D['endpoints'][1]), int(D['indices'][2]))
                      for D in triangulation]
        self.vertices = [Label(*v) for v in set(vertices)]
        
    def draw(self, shifts):
        glRasterPos3f(0.0, 0.0, 0.0)
        for M, L in shifts:
            disp = M*self.meridian + L*self.longitude
            glPushMatrix()
            glTranslatef(disp.real, disp.imag, 0.0)
            for labels in (self.segments, self.vertices):
                for label in labels:
                    label.draw()
            glPopMatrix()

cdef class HoroballScene:
    """
    A family of translations of a Horoball Group which fill the
    screen.  The horoballs are viewed by an observer sitting on one of
    the horoballs.  The variable which_cusp selects which cusp the
    viewer's horoball corresponds to.
    """
    cdef nbhd
    cdef meridian, longitude, offset, flipped
    cdef cusp_view, Ford, tri, pgram, labels, shifts
    cdef pgram_var, Ford_var, tri_var, horo_var, label_var
    cdef GLfloat Xangle, Yangle
    cdef GLuint ball_list_id, pgram_list_id, labels_list_id
    cdef GLuint tri_light_list_id, tri_dark_list_id
    cdef GLuint Ford_dark_list_id,  Ford_light_list_id
    cdef double cutoff
    cdef int which_cusp

    def __init__(self, nbhd, pgram_var, Ford_var, tri_var, horo_var, label_var,
                 flipped=False, cutoff=0.1, which_cusp=0):
        self.nbhd = nbhd
        self.which_cusp = which_cusp
        self.flipped = flipped
        self.tri_var = tri_var
        self.Ford_var = Ford_var
        self.pgram_var = pgram_var
        self.horo_var = horo_var
        self.label_var = label_var
        self.offset = 0.0j
        self.Xangle, self.Yangle = 0.0, 0.0
        self.pgram_list_id = base = glGenLists(7)
        self.ball_list_id = base + 1
        self.Ford_light_list_id = base + 2
        self.Ford_dark_list_id = base + 3
        self.tri_light_list_id = base + 4
        self.tri_dark_list_id = base + 5
        self.labels_list_id = base + 6
        self.set_cutoff(cutoff)
        self.pgram = Parallelogram()
        self.build_scene()

    def destroy(self):
        if self.cusp_view:
            self.cusp_view.delete_lists()
        if self.pgram_list_id > 0:
            glDeleteLists(self.pgram_list_id, 7)
            self.pgram_list_id = 0

    def set_cutoff(self, cutoff):
        self.cutoff = cutoff
        
    def flip(self, boolean_value):
        self.flipped = boolean_value

    def build_scene(self, which_cusp=None, full_list=True):
        if self.nbhd is None:
            self.cusp_view = self.Ford = self.tri = self.labels = None
            return
        if which_cusp == None:
            which_cusp = self.which_cusp
        else:
            self.which_cusp = which_cusp
        self.meridian, self.longitude = (
            complex(z) for z in self.nbhd.translations(self.which_cusp))
        self.cusp_view = HoroballGroup(
            self.nbhd.horoballs(self.cutoff, which_cusp, full_list),
            [self.nbhd.original_index(n) for n in range(self.nbhd.num_cusps())],
            self.meridian,
            self.longitude)
        self.Ford = FordEdgeSet(
                self.nbhd.Ford_domain(self.which_cusp),
                self.longitude, self.meridian)
        self.tri = TriangulationEdgeSet(
                self.nbhd.triangulation(self.which_cusp),
                self.longitude, self.meridian)
        self.labels = LabelSet(
                self.nbhd.triangulation(self.which_cusp),
                self.longitude, self.meridian)
        self.gl_compile()

    def build_shifts(self, R, T):
        self.shifts = []
        if self.cusp_view is None:
            return
        M = 1 + int(ceil(T/abs(self.meridian.imag)))
        N = 1 + int(ceil(R/self.longitude.real))
        for m in range(-M,M+1):
            shear = m*self.meridian.real/self.longitude.real
            left = int(floor(-shear-N))
            for n in range(left,left+2*N+1):
                self.shifts.append((m,n))

    def translate(self, z):
        """
        Translate modulo the cusp stabilizer.
        """
        if self.cusp_view is None:
            return
        if self.flipped:
            z = z.conjugate()
        z += self.offset
        z += 0.5*self.meridian.imag*1j
        z = z - (z.imag//self.meridian.imag)*self.meridian
        z -= 0.5*self.meridian.imag*1j
        z += 0.5*self.longitude.real
        z = z - (z.real//self.longitude.real)*self.longitude
        z -= 0.5*self.longitude
        self.offset = z
    
    cdef right_top(self):
        cdef GLfloat proj[16]
        glGetFloatv(GL_PROJECTION_MATRIX, proj)
        return (abs(<float>(1.0/proj[0])), abs(<float>(1.0/proj[5])))

    def gl_compile(self):
        self.pgram.build_display_list(self.pgram_list_id,
                                      self.longitude, self.meridian)
        right, top = self.right_top()
        R = right + 2.0 + 0.5*self.longitude.real
        T = top + 2.0 + 0.5*self.meridian.imag
        self.cusp_view.build_display_list(self.ball_list_id, R, T)
        self.build_shifts(R, T)
        self.Ford.build_display_list(self.Ford_light_list_id, self.shifts, dark=False)
        self.Ford.build_display_list(self.Ford_dark_list_id, self.shifts, dark=True)
        self.tri.build_display_list(self.tri_light_list_id, self.shifts, dark=False)
        self.tri.build_display_list(self.tri_dark_list_id, self.shifts, dark=True)
        self.labels.build_display_list(self.labels_list_id, self.shifts)

    def draw_segments(self, ford_height, pgram_height):
        with_horoballs = self.horo_var.get()
        glPushMatrix()
        glTranslatef(self.offset.real, self.offset.imag, ford_height)
        if self.tri_var.get():
            if with_horoballs:
                glCallList(self.tri_light_list_id)
            else:
                glCallList(self.tri_dark_list_id)
        if self.Ford_var.get():
            if with_horoballs:
                glCallList(self.Ford_light_list_id)
            else:
                glCallList(self.Ford_dark_list_id)
        glPopMatrix()
        if self.pgram_var.get():
            glPushMatrix()
            glTranslatef(0.0, 0.0, pgram_height)
            glCallList(self.pgram_list_id)
            glPopMatrix()

    def draw(self, *args):
        """
        The scene is drawn translated by self.offset, but the
        parallelogram stays fixed.
        """
        if self.nbhd is None:
            return
        glPushMatrix()
        if self.flipped:
            self.draw_segments(-2.0, -2.2)
            label_height = -2.4
        else:
            self.draw_segments(2.0, 2.2)
            label_height = 2.4
        if self.horo_var.get():
            glPushMatrix()
            glTranslatef(self.offset.real, self.offset.imag, 0.0)
            glCallList(self.ball_list_id)
            glPopMatrix()
        if self.label_var.get():
            glPushMatrix()
            glTranslatef(self.offset.real, self.offset.imag, label_height)
            glCallList(self.labels_list_id)
            glPopMatrix()
        glPopMatrix()

# Methods to translate and rotate our scene.

cdef glTranslateScene(s, x, y, mousex, mousey):
    cdef GLdouble X, Y
    cdef GLdouble mat[16]

    X, Y = s * (x - mousex), s * (mousey - y)
    glGetDoublev(GL_MODELVIEW_MATRIX, mat)
    glLoadIdentity()
    glTranslatef(X, Y, 0.0)
    glMultMatrixd(mat)

cdef glRotateScene(xcenter, ycenter, zcenter, Xangle, Yangle):
    cdef GLdouble mat[16]

    glGetDoublev(GL_MODELVIEW_MATRIX, mat)
    glLoadIdentity()
    glTranslatef(xcenter, ycenter, zcenter)
    glRotatef(Yangle, 1., 0., 0.)
    glRotatef(Xangle, 0., 1., 0.)
    glTranslatef(-xcenter, -ycenter, -zcenter)
    glMultMatrixd(mat)

class RawOpenGLWidget(Tk_.Widget, Tk_.Misc):
    """
    Widget without any sophisticated bindings
    by Tom Schwaller
    """

    def __init__(self, master, cnf={}, **kw):
        curr_platform = sys.platform
        cpu_width = platform.architecture()[0]
        if curr_platform[:5] == 'linux':
            curr_platform = 'linux2'
        if curr_platform[:5] == "linux" and cpu_width == '64bit':
            curr_platform += "-x86_64"
        if curr_platform == 'win32':
            windows_version = sys.getwindowsversion()
            if (windows_version.major, windows_version.minor) > (6,0):
                curr_platform += 'VC'
            if cpu_width == '64bit':
                curr_platform += '-x86_64'
        suffix = curr_platform + "-tk" + master.getvar("tk_version")
        Togl_path = os.path.abspath(os.path.join(togl.__path__[0], suffix))
        if not os.path.exists(Togl_path):
            raise RuntimeError('Togl directory "%s" missing.' % Togl_path)
        
        master.tk.call('lappend', 'auto_path', Togl_path)
        try:
            master.tk.call('package', 'require', 'Togl')
        except Tk_.TclError:
            raise RuntimeError('Tcl can not find Togl even though directory %s exists' % Togl_path)

        Tk_.Widget.__init__(self, master, 'togl', cnf, kw)
        self.root = master
        self.bind('<Map>', self.tkMap)
        self.bind('<Expose>', self.tkExpose)
        self.bind('<Configure>', self.tkExpose)

    def tkRedraw(self, *dummy):
        self.update_idletasks()
        self.tk.call(self._w, 'makecurrent')
        glPushMatrix()
        self.redraw()
        glPopMatrix()

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
                 mouse_scale=False,
                 fovy=30.0, near=1.0, far=100.0,
                 cnf={}, **kw):
        """
        Create an opengl widget.  Arrange for redraws when the window is
        exposed or when it changes size.
        """
        RawOpenGLWidget.__init__(*(self, master, cnf), **kw)
        self.help_text = help
        self.initialised = 0
        if sys.platform == 'darwin':
            self.config(cursor='hand')
        else:
            self.config(cursor='fleur')
        self.flipped = False

        # Current coordinates of the mouse.
        self.xmouse = self.ymouse = self.tmouse = self.delta_t = 0
        self.Xangle = self.Yangle = 0

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
        self.fovy = fovy

        # Position of clipping planes.
        self.near = near
        self.far = far

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
            self.bind('<B1-Motion>', self.tkTranslate)
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
            self.bind('<B1-Motion>', self.tkTranslate)
        if mouse_scale:
            self.bind('<Button-2>', self.tkRecordMouse)
            self.bind('<B2-Motion>', self.tkScale)
            self.bind('<KeyPress>', self.tkKeyPress)

    def help(self):
        """
        Help message for the widget.
        """
        InfoDialog(self, 'Viewer Help', self.help_text)

    def activate(self):
        """
        Cause this OpenGLWidget to be the current destination for
        drawing.  Does NOT make the widget be the focus of keyboard
        events; SnapPy OpenGL widgets do not accept keyboard events.
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
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        if redraw:
            self.tkRedraw()

    def tkHandlePick(self, event):
        """
        Handle a pick on the scene.
        """
        cdef GLdouble objX, objY, objZ
        cdef GLdouble model[16]
        cdef GLdouble proj[16]
        cdef GLint view[4]

        if hasattr(self, 'pick'):
            raise ValueError('Sorry, this implementation was removed')
            # See GLU/README.txt for how to restore.
            #
            # # here we need to use glu.UnProject
            # # Tk and X have their origin top left, 
            # # while OpenGLWidget has its origin bottom left.
            # # So we need to subtract y from the window height to get
            # # the proper pick position for OpenGLWidget
            # realy = self.winfo_height() - event.y
            # self.activate()
            # glGetDoublev(GL_MODELVIEW_MATRIX, model)
            # glGetDoublev(GL_PROJECTION_MATRIX, proj)
            # glGetIntegerv(GL_VIEWPORT, view)
            # gluUnProject(event.x, realy, 0., model, proj, view, &objX, &objY, &objZ)
            # p1 = (objX, objY, objZ)
            # gluUnProject(event.x, realy, 1., model, proj, view, &objX, &objY, &objZ)
            # p2 = (objX, objY, objZ)

            # if self.pick(self, p1, p2):
            #     # If the pick method returns true we redraw the scene.
            #     self.tkRedraw()

    def tkRecordMouse(self, event):
        """
        Record the current mouse position.
        """
        self.delta_t = event.time - self.tmouse
        self.xmouse, self.ymouse, self.tmouse = event.x, event.y, event.time

    def StartRotate(self, event):
        # Switch off any autospinning if it was happening
        self.autospin = 0
        self.tkRecordMouse(event)

    def tkScale(self, event):
        """
        Scale the scene.  Achieved by moving the eye position
        when using perspective.
        """
        scale = 1 - 0.01 * (event.y - self.ymouse)
        self.distance = self.distance * scale
        self.tkRedraw()
        self.tkRecordMouse(event)

    def zoom(self, x):
        t = float(x)/100.0
        self.distance = t*2.0 + (1-t)*10.0
        self.tkRedraw()

    def do_AutoSpin(self):
        self.activate()
        glRotateScene(self.xcenter, self.ycenter, self.zcenter,
                      self.Xangle, self.Yangle)
        self.tkRedraw()

        if self.autospin:
            self.after(10, self.do_AutoSpin)

    def tkAutoSpin(self, event):
        """
        Perform autospin of scene.
        """
        if self.autospin_allowed and 0 < self.delta_t < 100:
            self.autospin = 1
            self.after(10, self.do_AutoSpin)
        self.update_idletasks()

    def tkRotate(self, event):
        """
        Perform rotation of scene.
        """
        cdef GLfloat Xangle, Yangle
        self.activate()
        self.Xangle = 0.5 * (event.x - self.xmouse)
        self.Yangle = 0.5 * (event.y - self.ymouse)
        glRotateScene(self.xcenter, self.ycenter, self.zcenter,
                      self.Xangle, self.Yangle)
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

    def mouse_update(self, event):
        """
        Redraw the scene and save the mouse coordinates.
        """
        self.tkRedraw()
        self.tkRecordMouse(event)

    def tkRedraw(self, *dummy):
        """
        Cause the opengl widget to redraw itself.
        """
        if not self.initialised: return
        self.tk.call(self._w, 'makecurrent')
        glPushMatrix()                        # Protect our matrix
        self.update_idletasks()
        w = self.winfo_width()
        h = self.winfo_height()
        glViewport(0, 0, w, h)

        # Clear the background and depth buffer.
        glClearColor(self.r_back, self.g_back, self.b_back, 0.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        # build the projection matrix
        self.build_projection(w, h)

        # Call objects redraw method.
        self.redraw()
        glPopMatrix()                            # Restore the matrix
        self.tk.call(self._w, 'swapbuffers')

    def build_projection(self, width, height):
        cdef GLdouble xmax, yymax, near, far
        aspect = float(width)/float(height)
        near, far = self.near, self.far
        self.activate()
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        ymax = near * tan(self.fovy*pi/360.0)
        xmax = ymax * aspect
        glFrustum(-xmax, xmax, -ymax, ymax, near, far)
        glTranslatef(-self.xcenter, -self.ycenter, -(self.zcenter+self.distance))
        glMatrixMode(GL_MODELVIEW)

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
        # DEAL WITH THIS

class OpenGLOrthoWidget(OpenGLWidget):
    """
    A version of the widget that uses orthographic projection instead
    of perspective.
    """

    def build_projection(self, width, height, t=1.0):
        aspect = float(width)/float(height)
        top = self.fovy/2
        right = top*aspect
        self.activate()
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        if self.flipped:
            glEnable(GL_LIGHT1)
            glDisable(GL_LIGHT0)
            glOrtho(-right, right, top, -top, 3.0, -3.0)
        else:
            glEnable(GL_LIGHT0)
            glDisable(GL_LIGHT1)
            glOrtho(-right, right, -top, top, -3.0, 3.0)
        glMatrixMode(GL_MODELVIEW)

    def tkTranslate(self, event):
        """
        Perform translation of scene.
        """
        self.tkRedraw()
        self.tkRecordMouse(event)

    def redraw(self):
        pass
