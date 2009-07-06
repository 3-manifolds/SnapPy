#!/usr/bin/env python

import Tkinter
from Tkinter import * 
from OpenGL.GL import *
import OpenGL.GLU as GLU
from oglidget import Opengl
from colorsys import hls_to_rgb

import OpenGL, os, sys

def norm(vector):
  return sqrt(dot(vector, vector))

# This assumes there is only one cusp at the moment.
# ball_dicts is a list of lists of dicts, one list for each
# cusp and one dict for each fundamental horoball.
class HoroballViewer:

  def __init__(self, cusp_list, translation_list, root=None,
               title='Horoball Viewer'):
    self.title = title
    if root is None:
      root = Tkinter._default_root
    self.window = window = Toplevel(root)
    window.protocol("WM_DELETE_WINDOW", self.close)
    window.title(title)
    self.widget = widget = Opengl(master=self.window,
                                  width = 600,
                                  height = 600,
                                  double = 1,
                                  depth = 1,
                                  help = """
    XXX
""")
    widget.set_eyepoint(5.0)
    self.cusps = []
    for n in range(len(cusp_list)):
      self.cusps.append(HoroballGroup(cusp_list[n], translation_list[n]))
    widget.redraw = self.redraw
    widget.autospin_allowed = 0
    widget.set_background(.4, .4, .9)
    self.topframe = topframe = Frame(self.window, borderwidth=0,
                                     relief=FLAT, background='#f4f4f4')
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
    self.init_GL()
    self.init_matrix()
    self.set_lighting()

  # Subclasses may override this, e.g. if they use a help menu.
  def add_help(self):
    help = Button(self.topframe, text = 'Help', width = 4,
                  borderwidth=0, highlightthickness=0,
                  background="#f4f4f4", command = self.widget.help)
    help.grid(row=0, column=4, sticky=E, pady=3)
    self.topframe.columnconfigure(3, weight = 1)

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
    glFrontFace(GL_CCW);
    glMaterial(GL_FRONT, GL_AMBIENT, [1.0, 1.0, 1.0, 1.0])
    glMaterial(GL_FRONT, GL_DIFFUSE, [0.5, 0.5, 0.5, 1.0])
    glMaterial(GL_FRONT, GL_SPECULAR, [0.75, 0.75, 0.75, 1.0])
    glMaterial(GL_FRONT, GL_SHININESS, 100.0)
    
#    glMaterial(GL_BACK,  GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
#    glMaterial(GL_BACK,  GL_DIFFUSE, [0.5, 0.5, 0.5, 1.0])
#    glMaterial(GL_BACK,  GL_SPECULAR, [0.75, 0.75, 0.75, 1.0])
#    glMaterial(GL_BACK,  GL_SHININESS, 0.0)
    glEnable(GL_CULL_FACE)
    glCullFace(GL_BACK)

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
    self.init_matrix()  
    self.widget.set_eyepoint(10.0)
    self.zoom.set(50)
    self.widget.tkRedraw()

  def set_zoom(self, x):
    t = float(x)/100.0
    self.widget.distance = t*2.0 + (1-t)*20.0
    self.widget.tkRedraw()

  def redraw(self, widget):
    for cusp in self.cusps:
      cusp.draw(widget)

class HoroballGroup:
  """
  A fundamental set of horoballs for a single cusp.
  """
  def __init__(self, dicts, translations):
    self.display_list = glGenLists(1)
    self.domain_list = glGenLists(1)
    self.meridian, self.longitude = translations
    self.build_spheres(dicts)
    glNewList(self.display_list, GL_COMPILE)
    for m in range(-1,2):
      for n in range(-1,2):
        translation = m*self.meridian + n*self.longitude
        glPushMatrix()
        glTranslate(translation.real, translation.imag, 0.0)
        glCallList(self.domain_list)
        glPopMatrix()
    glEndList()    

  def draw(self, widget):
    glCallList(self.display_list)

  def build_spheres(self, dicts):
    self.sphere_quadric = GLU.gluNewQuadric()
    GLU.gluQuadricDrawStyle(self.sphere_quadric, GLU.GLU_FILL)
    GLU.gluQuadricNormals(self.sphere_quadric, GLU.GLU_SMOOTH)
    glNewList(self.domain_list, GL_COMPILE) 
    for D in dicts:
      center = [D['center'].real, D['center'].imag, D['radius']]
      radius = D['radius']
      # the color should depend on the cusp index
      glPushMatrix()
      glTranslate(*center)
      glColor4f(1.0, 0.2, 0.2, .8)
      GLU.gluSphere(self.sphere_quadric, radius, 20, 20)
      glPopMatrix()
    glEndList()

__doc__ = """
   The horoviewer module exports the HoroballViewer class, which is
   a Tkinter / OpenGL window for viewing cusp neighborhoods.
   """

__all__ = ['HoroballViewer']

# data for testing
test_cusps =[ [
{'index': 0, 'radius': 0.10495524653309458, 'center': (-0.25834708978942406+1.4921444888317912j)},
{'index': 0, 'radius': 0.10495524653309474, 'center': (-1.9740640711918203+2.2591328216964328j)},
{'index': 0, 'radius': 0.10495524653309476, 'center': (-0.85785849070119979+0.38349416643232093j)},
{'index': 0, 'radius': 0.10495524653309482, 'center': (-0.39768176782278597-0.85240383024834776j)},
{'index': 0, 'radius': 0.1142737629103735, 'center': (-1.6448178786236118+2.0647270360966177j)},
{'index': 0, 'radius': 0.11427376291037381, 'center': (-0.58759328235763275+1.6865502744316054j)},
{'index': 0, 'radius': 0.11427376291037425, 'center': (-0.72692796039099428-0.65799804464853373j)},
{'index': 0, 'radius': 0.11427376291037446, 'center': (-0.5286122981329916+0.18908838083250662j)},
{'index': 0, 'radius': 0.16407530192719913, 'center': (-1.6831397382358244+1.7935639049681977j)},
{'index': 0, 'radius': 0.16407530192719913, 'center': (-1.8048116812694035+1.4888037417440001j)},
{'index': 0, 'radius': 0.16407530192719999, 'center': (-0.56693415774520339-0.082074750295914087j)},
{'index': 0, 'radius': 0.16407530192719999, 'center': (-0.68860610077878237-0.38683491352011279j)},
{'index': 0, 'radius': 0.20193197944608232, 'center': (-0.28353891189056118-0.47924212868872834j)},
{'index': 0, 'radius': 0.20193197944608232, 'center': (-2.0882069271240447+1.8859711201368135j)},
{'index': 0, 'radius': 0.2019319794460824, 'center': (0.28353891189056052+0.47924212868872862j)},
{'index': 0, 'radius': 0.2019319794460824, 'center': (-1.3997444923811837+1.3963965265753839j)},
{'index': 0, 'radius': 0.27835650978783572, 'center': (-0.76862048805880945+1.3219220338932829j)},
{'index': 0, 'radius': 0.27835650978783572, 'center': (-0.34758509243181374+0.55371662137082966j)},
{'index': 0, 'radius': 0.27835650978783577, 'center': (-0.90795516609217164-1.0226262851868568j)},
{'index': 0, 'radius': 0.27835650978783577, 'center': (-2.7193309314464176+1.9604456128189169j)},
{'index': 0, 'radius': 0.38387596322871925, 'center': 0j},
{'index': 0, 'radius': 0.38387596322871925, 'center': (-2.3717458390146056+1.4067289914480869j)},
{'index': 0, 'radius': 0.5, 'center': (0.069667339016679999+1.1722741595400699j)},
{'index': 0, 'radius': 0.5, 'center': (-2.3020784999979229+2.5790031509881559j)}] ]

test_translations = [( (1.2555402585239854+0.46890966381602794j),
                       (12.276733229173129+0j) )]

if __name__ == '__main__':
    HV = HoroballViewer(test_cusps, test_translations)
    HV.window.mainloop()


