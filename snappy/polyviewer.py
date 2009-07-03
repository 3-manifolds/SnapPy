from snappy.CyOpenGL import *
import Tkinter

class PolyhedronViewer:
    """
    Window for viewing a hyperbolic polyhedron, either in the Poincare
    or Klein model.
    """

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
  the polyhedron if you zoom far enough.
""")
        widget.set_eyepoint(5.0)
        self.model_var=StringVar(value='Klein')
        self.sphere_var=IntVar(value=1)
        self.GL = GL_context()
        self.polyhedron = HyperbolicPolyhedron(facedicts,
                                               self.model_var,
                                               self.sphere_var)
        widget.redraw = self.polyhedron.draw
        widget.autospin_allowed = 1
        widget.set_background(.2, .2, .2)
        self.topframe = topframe = Frame(self.window, borderwidth=0,
                                         relief=FLAT, background='#f4f4f4')
        self.klein = Radiobutton(topframe, text='Klein', value='Klein',
                                 variable = self.model_var,
                                 command = self.new_model,
                                 background='#f4f4f4')
        self.poincare = Radiobutton(topframe, text=u'Poincar\u00e9', value='Poincare',
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
