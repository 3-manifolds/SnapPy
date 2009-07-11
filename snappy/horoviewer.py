#!/usr/bin/env python

import Tkinter
from Tkinter import * 
from snappy.CyOpenGL import *
from colorsys import hls_to_rgb

import os, sys

class HoroballViewer: 

  def __init__(self, cusp_list, translation_list, Ford_segments, triangulation,
               which_cusp=0, root=None, title='Horoball Viewer'):
    self.title = title
    if root is None:
      root = Tkinter._default_root
    self.window = window = Toplevel(root)
    window.protocol("WM_DELETE_WINDOW", self.close)
    window.title(title)
    self.widget = widget = OpenGLWidget(master=self.window,
                                        width = 600,
                                        height = 600,
                                        double = 1,
                                        depth = 1,
                                        help = """
    XXX
""")
    widget.set_eyepoint(5.0)
    self.pgram_var = pgram_var = Tk_.IntVar(value=1)
    self.Ford_var = Ford_var = Tk_.IntVar(value=1)
    self.tri_var = tri_var = Tk_.IntVar(value=1)
    self.GL = GL_context()
    self.GLU = GLU_context()
    self.scene = HoroballScene(cusp_list, translation_list, Ford_segments,
                               triangulation, pgram_var, Ford_var, tri_var,
                               which_cusp)
    widget.redraw = self.scene.draw
    widget.autospin_allowed = 0
    widget.set_background(.5, .5, .5)
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

  def reset(self):
#    self.init_matrix()  
    self.widget.set_eyepoint(10.0)
    self.zoom.set(50)
    self.widget.tkRedraw()

  def set_zoom(self, x):
    t = float(x)/100.0
    self.widget.distance = t*2.0 + (1-t)*20.0
    self.widget.tkRedraw()

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


