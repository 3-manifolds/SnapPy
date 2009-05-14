#!/usr/bin/env python

from Tkinter import * 
from OpenGL.GL import *
from OpenGL.GLUT import glutSolidSphere
from oglidget import Opengl
from colorsys import hls_to_rgb

import OpenGL, os, sys

def norm(vector):
  return sqrt(dot(vector, vector))

# This assumes there is only one cusp at the moment.
# ball_dicts is a list of lists of dicts, one list for each
# cusp and one dict for each fundamental horoball.
class HoroballViewer:

  def __init__(self, cusp_list, translation_list, title='Horoball Viewer'):
    self.window = Tk()
    self.window.title(title)
    self.widget = Opengl(master=self.window,
                         width = 600,
                         height = 600,
                         double = 1,
                         depth = 1,
                         help = """
    XXX
""")
    self.widget.set_eyepoint(5.0)
    self.cusps = []
    for n in range(len(cusp_list)):
      self.cusps.append(HoroballGroup(cusp_list[n], translation_list[n]))
    self.widget.redraw = self.redraw
    self.widget.autospin_allowed = 0
    self.widget.set_background(.4, .4, .9)
    frame       = Frame(self.window, bd = 3, relief = RIDGE)
    quit  = Button(frame, text = 'Quit', pady = 1, width = 4, 
                     command = self.window.destroy)
    reset = Button(frame, text = 'Reset', pady = 1, width = 4,
                     command = self.reset)
    help  = Button(frame, text = 'Help', pady = 1, width = 4,
                     command = self.widget.help)
    self.zoom = Scale(self.window, showvalue=0, from_=100, to=0,
           command = self.set_zoom)
    self.zoom.set(50)

    quit.grid(row=0, column=0)
    reset.grid(row=0, column=1, sticky=W)
    help.grid(row=0, column=2, sticky=W)
    frame.columnconfigure(2, weight = 1)
    frame.pack(side = TOP, fill = X)
    self.widget.pack(side = LEFT, expand = YES, fill = BOTH)
    self.zoom.pack(side = RIGHT, fill = Y)
    self.widget.extra_help = 'HELP'
    self.init_GL()
    self.init_matrix()
    self.set_lighting()

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
    glNewList(self.domain_list, GL_COMPILE) 
    for dict in dicts:
      center = [dict['center'].real, dict['center'].imag, dict['radius']]
      radius = dict['radius']
      # the color should depend on the cusp index
      glPushMatrix()
      glTranslate(*center)
      glColor4f(1.0, 0.2, 0.2, .8)
      glutSolidSphere(radius, 20, 20)
      glPopMatrix()
    glEndList()

__doc__ = """
   The polyviewer module exports the HoroballViewer class, which is
   a Tkinter / OpenGL window for viewing cusp neighborhoods.
   """

__all__ = ['HoroballViewer']

# data for testing
testdata =[{'index': 0, 'radius': 0.017187631363731467, 'center': (0.7629506355167095+0.042188120412997687j)},
 {'index': 0, 'radius': 0.017187631363731467, 'center': (-0.76295063551670994-0.33026662618317804j)},
 {'index': 0, 'radius': 0.017187631363731474, 'center': (0.92804428382437809-0.042188120412999963j)},
 {'index': 0, 'radius': 0.017187631363731474, 'center': (0.19928566240301399+0.042188120412998437j)},
 {'index': 0, 'radius': 0.017187631363731481, 'center': (-0.36437931071068258-0.24589038535718008j)},
 {'index': 0, 'radius': 0.017187631363731481, 'center': (0.36437931071068219-0.042188120412999491j)},
 {'index': 0, 'radius': 0.017187631363731491, 'center': (1.0447831220735575+0.10185113247208966j)},
 {'index': 0, 'radius': 0.017187631363731491, 'center': (0.6462117972675302-0.10185113247209188j)},
 {'index': 0, 'radius': 0.017187631363731491, 'center': (-0.082546824153834658-0.18622737329808842j)},
 {'index': 0, 'radius': 0.017187631363731491, 'center': (-0.48111814895986199-0.38992963824226967j)},
 {'index': 0, 'radius': 0.017187631363731498, 'center': (-1.0447831220735582-0.38992963824227117j)},
 {'index': 0, 'radius': 0.017187631363731498, 'center': (0.48111814895986205+0.1018511324720901j)},
 {'index': 0, 'radius': 0.017187631363731505, 'center': (0.082546824153834242-0.10185113247209115j)},
 {'index': 0, 'radius': 0.017187631363731505, 'center': (-0.64621179726753075-0.18622737329808811j)},
 {'index': 0, 'radius': 0.017187631363731512, 'center': (-0.92804428382437898-0.24589038535718094j)},
 {'index': 0, 'radius': 0.017187631363731512, 'center': (-0.19928566240301396-0.042188120412998451j)}]


if __name__ == '__main__':
    HV = HoroballViewer(testdata)
    HV.window.mainloop()


