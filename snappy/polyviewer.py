from __future__ import unicode_literals
# -*- coding: utf-8 -*-
from snappy.CyOpenGL import *

try:
    import Tkinter as Tk_
except ImportError: # Python 3
    import tkinter as Tk_

class PolyhedronViewer:
    """
    Window for viewing a hyperbolic polyhedron, either in the Poincare
    or Klein model.
    """

    def __init__(self, facedicts, root=None, title='Polyhedron Viewer'):
        self.title=title
        if root is None:
            root = Tk_._default_root
        self.window = window = Tk_.Toplevel(master=root, class_='snappy')
        window.withdraw()
        window.title(title)
        window.protocol("WM_DELETE_WINDOW", self.close)
        self.topframe = topframe = Tk_.Frame(window, borderwidth=0,
                                             relief=Tk_.FLAT, background='#f4f4f4')
        self.bottomframe = bottomframe = Tk_.Frame(window, borderwidth=0,
                                             relief=Tk_.FLAT)
        self.widget = widget = OpenGLWidget(master=bottomframe,
                                            width = 600,
                                            height = 600,
                                            double = 1,
                                            depth = 1,
                                            help = """
  Use mouse button 1 to rotate the polyhedron.
  Releasing the button while moving will "throw"
  the polyhedron and make it keep spinning.

  The slider controls zooming.  You can see inside
  the polyhedron if you zoom far enough.
""")
        widget.set_eyepoint(5.0)
        self.model_var=Tk_.StringVar(value='Klein')
        self.sphere_var=Tk_.IntVar(value=1)
        self.GL = GL_context()
        self.polyhedron = HyperbolicPolyhedron(facedicts,
                                               self.model_var,
                                               self.sphere_var)
        widget.redraw = self.polyhedron.draw
        widget.autospin_allowed = 1
        widget.set_background(.2, .2, .2)

        self.klein = Tk_.Radiobutton(topframe, text='Klein',
                                     value='Klein',
                                     variable = self.model_var,
                                     command = self.new_model,
                                     background='#f4f4f4')
        self.poincare = Tk_.Radiobutton(topframe, text='Poincaré',
                                        value='Poincare',
                                        variable = self.model_var,
                                        command = self.new_model,
                                        background='#f4f4f4')
        self.sphere = Tk_.Checkbutton(topframe, text='',
                                      variable = self.sphere_var,
                                      command = self.new_model,
                                      borderwidth=0, background='#f4f4f4')
        self.spherelabel = Tk_.Text(topframe, height=1, width=3,
                                    relief=Tk_.FLAT, font='Helvetica 14 bold',
                                    borderwidth=0, highlightthickness=0,
                                    background='#f4f4f4')
        self.spherelabel.tag_config("sub", offset=-4)
        self.spherelabel.insert(Tk_.END, 'S')
        self.spherelabel.insert(Tk_.END, '∞', 'sub')
        self.spherelabel.config(state=Tk_.DISABLED)

        self.klein.grid(row=0, column=0, sticky=Tk_.W, padx=20)
        self.poincare.grid(row=0, column=1, sticky=Tk_.W, padx=20)
        self.sphere.grid(row=0, column=2, sticky=Tk_.W, padx=0)
        self.spherelabel.grid(row=0, column=3, sticky=Tk_.W)
        self.add_help()
        topframe.pack(side=Tk_.TOP, fill=Tk_.X)
        self.zoomframe = zoomframe = Tk_.Frame(bottomframe, borderwidth=0, relief=Tk_.FLAT)
        self.zoom = zoom = Tk_.Scale(zoomframe, showvalue=0, from_=100, to=0,
                                     command = self.set_zoom, width=11,
                                     troughcolor='#f4f4f4', borderwidth=1,
                                     relief=Tk_.SUNKEN)
        zoom.set(50)
        spacer = Tk_.Frame(zoomframe, height=14, borderwidth=0, relief=Tk_.FLAT)
        zoom.pack(side=Tk_.TOP, expand=Tk_.YES, fill=Tk_.Y)
        spacer.pack()
        bottomframe.columnconfigure(0, weight=1)
        widget.grid(row=0, column=0, sticky=Tk_.EW)
        zoomframe.grid(row=0, column=1, sticky=Tk_.NS)
        bottomframe.pack(side=Tk_.TOP, expand=Tk_.YES, fill=Tk_.BOTH)
        self.build_menus()
        window.deiconify()
        window.update() # Seems to avoid a race condition with togl
        self.bottomframe.bind('<Configure>', self.togl_handle_resize)

  # Subclasses may override this, e.g. if there is a help menu already.
    def add_help(self):
        help = Tk_.Button(self.topframe, text = 'Help', width = 4,
                          borderwidth=0, highlightthickness=0,
                          background="#f4f4f4", command = self.widget.help)
        help.grid(row=0, column=4, sticky=Tk_.E, pady=3)
        self.topframe.columnconfigure(3, weight = 1)
        #self.widget.extra_help = 'HELP'

  # Subclasses may override this to provide menus.
    def build_menus(self):
        pass

    def close(self):
        self.polyhedron.destroy()
        self.window.destroy()
        
    def reset(self):
        self.widget.autospin = 0
#        self.init_matrix()  
        self.widget.set_eyepoint(5.0)
        self.zoom.set(50)
        self.widget.tkRedraw()

    def set_zoom(self, x):
        t = float(x)/100.0
        self.widget.distance = t*1.0 + (1-t)*8.0
        self.widget.tkRedraw()

    def new_model(self):
        self.widget.tkRedraw()

    def togl_handle_resize(self, event):
        self.widget.config(height=self.bottomframe.winfo_height())
        self.widget.redraw()

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
