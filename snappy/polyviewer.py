#!/usr/bin/env python

# If this is too slow, we'll think about using Cython.

import Tkinter, OpenGL, os, sys
from Tkinter import * 
from colorsys import hls_to_rgb
from OpenGL.GL import *
from OpenGL.GLUT import glutInit, glutWireSphere
from oglidget import Opengl
from math import sqrt

class vector3:
    """
    A simple real 3-dimensional vector which supports addition,
    subtraction and right multiplication or division by scalars.
    Attributes include its norm and the square of its norm.
    """
    def __init__(self, triple):
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

class PolyhedronViewer:

  def __init__(self, facedicts, root=None, title=u'Polyhedron Viewer'):
    self.title=title
    if root is None:
      root = Tkinter._default_root
    self.window = window = Toplevel(root)
    window.title(title)
    window.protocol("WM_DELETE_WINDOW", self.close)
    self.widget = widget = Opengl(master=self.window,
                                  width = 600,
                                  height = 600,
                                  double = 1,
                                  depth = 1,
                                  help = """
  Use mouse button 1 to rotate the polyhedron.
  Releasing the button while moving will "throw"
  the polyhedron and make it keep spinning.

  The slider controls zooming.  You can see inside
  the polyhedron if you soom far enough.
""")
    widget.set_eyepoint(5.0)
    self.model_var=StringVar(value='Klein')
    self.sphere_var=IntVar(value=1)
    self.init_GL()
    self.init_matrix()
    self.set_lighting()
    if sys.platform != 'darwin':
        glutInit()
    self.polyhedron = HyperbolicPolyhedron(facedicts,
                                           self.model_var,
                                           self.sphere_var)
    widget.redraw = self.polyhedron.draw
    widget.autospin_allowed = 1
    widget.set_background(.4, .4, .9)
    self.topframe = topframe = Frame(self.window, borderwidth=0,
                                     relief=FLAT, background='#f4f4f4')
    self.klein = Radiobutton(topframe, text='Klein', value='Klein',
                             variable = self.model_var,
                             command = self.new_model,
                             background='#f4f4f4')
    self.poincare = Radiobutton(topframe, text='Poincare', value='Poincare',
                                variable = self.model_var,
                                command = self.new_model,
                                background='#f4f4f4')
    self.sphere = Checkbutton(topframe, text='',
                              variable = self.sphere_var,
                              command = self.new_model,
                              borderwidth=0, background='#f4f4f4')
    self.spherelabel = Text(topframe, height=1, width=3,
                            relief=FLAT, font='Helvetica 14 bold',
                            borderwidth=0, highlightthickness=0,
                            background='#f4f4f4')
    self.spherelabel.tag_config("sub", offset=-4)
    self.spherelabel.insert(END, 'S')
    self.spherelabel.insert(END, u'\u221e', "sub")
    self.spherelabel.config(state=DISABLED)

    self.klein.grid(row=0, column=0, sticky=W, padx=20)
    self.poincare.grid(row=0, column=1, sticky=W, padx=20)
    self.sphere.grid(row=0, column=2, sticky=W, padx=0)
    self.spherelabel.grid(row=0, column=3, sticky=W)
    self.add_help()
    topframe.pack(side=TOP, fill=X)
    widget.pack(side=LEFT, expand=YES, fill=BOTH)
    zoomframe = Frame(self.window, borderwidth=0, relief=FLAT)
    self.zoom = zoom = Scale(zoomframe, showvalue=0, from_=100, to=0,
                             command = self.set_zoom, width=11,
                             troughcolor='#f4f4f4', borderwidth=1,
                             relief=SUNKEN)
    zoom.set(50)
    spacer = Frame(zoomframe, height=14, borderwidth=0, relief=FLAT)
    zoom.pack(side=TOP, expand=YES, fill=Y)
    spacer.pack()
    zoomframe.pack(side=RIGHT, expand=YES, fill=Y)
    self.build_menus()


  # Subclasses may override this, e.g. if there is a help menu already.
  def add_help(self):
    help = Button(self.topframe, text = 'Help', width = 4,
                  borderwidth=0, highlightthickness=0,
                  background="#f4f4f4", command = self.widget.help)
    help.grid(row=0, column=4, sticky=E, pady=3)
    self.topframe.columnconfigure(3, weight = 1)
    #self.widget.extra_help = 'HELP'

  # Subclasses may override this to provide menus.
  def build_menus(self):
    pass

  def close(self):
      self.window.destroy()

  def init_GL(self):
    glEnable(GL_COLOR_MATERIAL)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_ALPHA_TEST)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
    glShadeModel(GL_SMOOTH)
    glEnable(GL_LIGHTING)
    glMatrixMode(GL_MODELVIEW);

  def set_lighting(self):
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, [0.5, 0.5, 0.5, 1.0])
    glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, 1.0)
    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, TRUE)
    glEnable(GL_LIGHT0)
    glLightfv(GL_LIGHT0, GL_POSITION, [2.0, 2.0, 2.0, 1.0] )
    glLightfv(GL_LIGHT0, GL_DIFFUSE, [0.5, 0.5, 0.5, 1.0] )
    glLightfv(GL_LIGHT0, GL_SPECULAR, [0.5, 0.5, 0.5, 1.0] )

  def init_matrix(self):
    glLoadIdentity()
#    glRotatef(30, 0.0, 0.0, -1.0)
#    glRotatef(45, -sqrt(3.0)*0.5, -0.5, 0.0)

  def reset(self):
    self.widget.autospin = 0
    self.init_matrix()  
    self.widget.set_eyepoint(5.0)
    self.zoom.set(50)
    self.widget.tkRedraw()

  def set_zoom(self, x):
    t = float(x)/100.0
    self.widget.distance = t*1.0 + (1-t)*8.0
    self.widget.tkRedraw()

  def new_model(self):
    self.widget.tkRedraw()

class PoincareTriangle:

  def __init__(self, vertices, center):
    self.vertices = vertices
    self.center = center

  def render_vertex(self, vertex):
    scale = 1 + sqrt(max(0, 1 - vertex.norm_squared))
    V = vertex/scale
    N = self.center - V
    N = N/N.norm
    glNormal3f(N.x, N.y, N.z)
    glVertex3f(V.x, V.y, V.z)

  def render(self, depth=5):
    C, V1, V2 = self.vertices
    M = (V1 + V2)/2
    CV1 = V1 - C
    MV1 = V1 - M
    CV2 = V2 - C
    MV2 = V2 - M
    step = 1.0/depth
    glBegin(GL_TRIANGLE_STRIP)
    for n in range(depth):
      self.render_vertex(M + MV1*n*step)
      self.render_vertex(C + CV1*n*step)
    self.render_vertex(V1)
    glEnd()
    glBegin(GL_TRIANGLE_STRIP)
    for n in range(depth):
      self.render_vertex(C + CV2*n*step)
      self.render_vertex(M + MV2*n*step)
    self.render_vertex(V2)
    glEnd()

class Face:
   """
   A face of a hyperbolic polyhedron.  Instantiate as
   Face(vertices=[...], closest=[x,y,z], distance=d, hue=h)
   vertices: list of vertex coordinates (in the Klein model)
   distance: distance from the origin to the plane
   closest: coordinates of the point nearest the origin on the plane
            containing the face
   hue: (float) hue to use in coloring the face.
   """
   def __init__(self, vertices, distance, closest, hue):
     self.vertices = [vector3(v) for v in vertices]
     self.distance = distance
     self.closest = vector3(closest)
     K = self.closest/self.closest.norm
     self.klein_normal = (K.x, K.y, K.z)
     self.center = self.closest/self.closest.norm_squared
     self.hue = hue

   def triangulate(self):
     Vlist = self.vertices
     zero = vector3((0,0,0))
     N = len(Vlist)
     triangle_list = []
     centroid = sum(Vlist, zero)/N
     for i in range(0,N):
       vertices = [centroid, Vlist[i-1],Vlist[i]]
       triangle_list.append(PoincareTriangle(vertices, self.center))
     return triangle_list

   def klein_render(self):
     r, g, b = hls_to_rgb(self.hue,0.5, 1.0) 
     glColor4f(r, g, b, 1.0)
     glBegin(GL_POLYGON)
     glNormal3f(*self.klein_normal)
     for V in self.vertices:
       glVertex3f(V.x, V.y, V.z)
     glEnd()

   def poincare_render(self):
     r, g, b = hls_to_rgb(self.hue, 0.5, 1.0) 
     glColor4f(r, g, b, 1.0)
     triangles = self.triangulate()
     for triangle in triangles:
       triangle.render()

class HyperbolicPolyhedron:
   """
   A hyperbolic polyhedron for display in OpenGL, either in the
   Klein model or the Poincare model.  Includes a representation
   of the sphere at infinity.
   """

   def __init__(self, facedicts, model_var, sphere_var):
     self.model = model_var
     self.sphere = sphere_var
     self.sphere_list = glGenLists(1)
     self.klein_list = glGenLists(1)
     self.poincare_list = glGenLists(1)
     self.faces = [Face(**dict) for dict in facedicts]
     self.build_sphere(self.sphere_list)
     self.build_klein_poly(self.klein_list)
     self.build_poincare_poly(self.poincare_list)

   def draw(self, widget):
     model = self.model.get()
     if model == 'Klein':
       glCallList(self.klein_list)
     elif model == 'Poincare':
       glCallList(self.poincare_list)
     if self.sphere.get():
       glPushMatrix()
       glLoadIdentity()
       glCallList(self.sphere_list)
       glPopMatrix()

   def set_poly_material(self):
     glMaterial(GL_FRONT, GL_AMBIENT, [1.0, 1.0, 1.0, 1.0])
     glMaterial(GL_FRONT, GL_DIFFUSE, [0.5, 0.5, 0.5, 1.0])
     glMaterial(GL_FRONT, GL_SPECULAR, [0.5, 0.5, 0.5, 1.0])
     glMaterial(GL_FRONT, GL_SHININESS, 50.0)
 
     glMaterial(GL_BACK,  GL_AMBIENT, [1.0, 1.0, 1.0, 1.0])
     glMaterial(GL_BACK,  GL_DIFFUSE, [0.5, 0.5, 0.5, 1.0])
     glMaterial(GL_BACK,  GL_SPECULAR, [0.5, 0.5, 0.5, 1.0])
     glMaterial(GL_BACK,  GL_SHININESS, 50.0)

   def build_sphere(self, list):
     glNewList(list, GL_COMPILE) 
     glRotatef(90, 1.0, 0.0, 0.0)
#     glEnable(GL_CULL_FACE);
#     glCullFace(GL_BACK);
     glFrontFace(GL_CCW);
     glMaterial(GL_FRONT, GL_AMBIENT, [1.0, 1.0, 1.0, 1.0])
     glMaterial(GL_FRONT, GL_DIFFUSE, [0.5, 0.5, 0.5, 1.0])
     glMaterial(GL_FRONT, GL_SPECULAR, [0.5, 0.5, 0.5, 1.0])
     glMaterial(GL_FRONT, GL_SHININESS, 50.0)
 
     glMaterial(GL_BACK,  GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
     glMaterial(GL_BACK,  GL_DIFFUSE, [0.5, 0.5, 0.5, 1.0])
     glMaterial(GL_BACK,  GL_SPECULAR, [0.5, 0.5, 0.5, 1.0])
     glMaterial(GL_BACK,  GL_SHININESS, 50.0)
     glColor4f(0.8, 0.8, 1.0, .2)
     glutWireSphere(1.0, 50, 50)
     glEndList()

   def build_klein_poly(self, list):
     glNewList(list, GL_COMPILE) 
     glDisable(GL_CULL_FACE);
#     glCullFace(GL_BACK);
#     glFrontFace(GL_CCW);
     self.set_poly_material()
     for face in self.faces:
       face.klein_render()
     glEndList()

   def build_poincare_poly(self, list):
     glNewList(list, GL_COMPILE) 
     glDisable(GL_CULL_FACE);
     self.set_poly_material()
     for face in self.faces:
       face.poincare_render()
     glEndList()

__doc__ = """
   The polyviewer module exports the PolyhedronViewer class, which is
   a Tkinter / OpenGL window for viewing Dirichlet Domains in either
   the Klein model or the Poincare model.
   """

__all__ = ['PolyhedronViewer']

# data for testing
testpoly = [{'distance': 0.57940518021497345,
 'vertices': [(0.34641016151377546, -0.34641016151377546, 0.34641016151377546), (0.57735026918962595, -0.57735026918962595, -0.57735026918962562), (0.57735026918962573, 0.57735026918962573, 0.57735026918962573)],
 'closest': [0.4723774929733302, -0.15745916432444337, 0.15745916432444337],
 'hue': 0.0},
 {'distance': 0.57940518021497345,
 'vertices': [(-0.57735026918962529, -0.57735026918962573, 0.57735026918962562), (-0.34641016151377557, -0.34641016151377557, -0.34641016151377518), (0.57735026918962595, -0.57735026918962595, -0.57735026918962562)],
 'closest': [-0.15745916432444337, -0.4723774929733302, -0.15745916432444337], 'hue': 0.5},
 {'distance': 0.57940518021497345, 'vertices': [(-0.34641016151377546, 0.34641016151377541, 0.34641016151377541), (0.57735026918962573, 0.57735026918962573, 0.57735026918962573), (-0.57735026918962573, 0.57735026918962573, -0.57735026918962573)],
 'closest': [-0.15745916432444337, 0.4723774929733302, 0.15745916432444337],
 'hue': 0.25}, {'distance': 0.57940518021497345,
 'vertices': [(0.57735026918962595, -0.57735026918962595, -0.57735026918962562), (-0.34641016151377557, -0.34641016151377557, -0.34641016151377518), (-0.57735026918962573, 0.57735026918962573, -0.57735026918962573)],
 'closest': [-0.15745916432444337, -0.15745916432444337, -0.4723774929733302],
 'hue': 0.25}, {'distance': 0.57940518021497345,
 'vertices': [(0.57735026918962573, 0.57735026918962573, 0.57735026918962573), (-0.57735026918962529, -0.57735026918962573, 0.57735026918962562), (0.34641016151377546, -0.34641016151377546, 0.34641016151377546)],
 'closest': [0.15745916432444337, -0.15745916432444337, 0.4723774929733302],
 'hue': 0.75}, {'distance': 0.57940518021497345,
 'vertices': [(-0.57735026918962573, 0.57735026918962573, -0.57735026918962573), (-0.34641016151377557, -0.34641016151377557, -0.34641016151377518), (-0.57735026918962529, -0.57735026918962573, 0.57735026918962562)],
 'closest': [-0.4723774929733302, -0.15745916432444337, -0.15745916432444337],
 'hue': 0.75}, {'distance': 0.57940518021497345,
 'vertices': [(0.57735026918962595, -0.57735026918962595, -0.57735026918962562), (0.34641016151377568, 0.34641016151377546, -0.34641016151377535), (0.57735026918962573, 0.57735026918962573, 0.57735026918962573)],
 'closest': [0.4723774929733302, 0.15745916432444337, -0.15745916432444337],
 'hue': 0.125}, {'distance': 0.57940518021497345,
 'vertices': [(-0.57735026918962529, -0.57735026918962573, 0.57735026918962562), (0.57735026918962573, 0.57735026918962573, 0.57735026918962573), (-0.34641016151377546, 0.34641016151377541, 0.34641016151377541)],
 'closest': [-0.15745916432444337, 0.15745916432444337, 0.4723774929733302],
 'hue': 0.125}, {'distance': 0.57940518021497345,
 'vertices': [(0.34641016151377546, -0.34641016151377546, 0.34641016151377546), (-0.57735026918962529, -0.57735026918962573, 0.57735026918962562), (0.57735026918962595, -0.57735026918962595, -0.57735026918962562)],
 'closest': [0.15745916432444337, -0.4723774929733302, 0.15745916432444337],
 'hue': 0.625}, {'distance': 0.57940518021497345,
 'vertices': [(0.57735026918962595, -0.57735026918962595, -0.57735026918962562), (-0.57735026918962573, 0.57735026918962573, -0.57735026918962573), (0.34641016151377568, 0.34641016151377546, -0.34641016151377535)],
 'closest': [0.15745916432444337, 0.15745916432444337, -0.4723774929733302],
 'hue': 0.625}, {'distance': 0.57940518021497345,
 'vertices': [(-0.57735026918962529, -0.57735026918962573, 0.57735026918962562), (-0.34641016151377546, 0.34641016151377541, 0.34641016151377541), (-0.57735026918962573, 0.57735026918962573, -0.57735026918962573)],
 'closest': [-0.4723774929733302, 0.15745916432444337, 0.15745916432444337],
 'hue': 0.0}, {'distance': 0.57940518021497345,
 'vertices': [(-0.57735026918962573, 0.57735026918962573, -0.57735026918962573), (0.57735026918962573, 0.57735026918962573, 0.57735026918962573), (0.34641016151377568, 0.34641016151377546, -0.34641016151377535)],
 'closest': [0.15745916432444337, 0.4723774929733302, -0.15745916432444337],
 'hue': 0.5}]

if __name__ == '__main__':
    PV = PolyhedronViewer(testpoly)
    PV.window.mainloop()


