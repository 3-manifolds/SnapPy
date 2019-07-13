# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from .CyOpenGL import *
from .export_stl import stl
from .theme import SnapPyStyle
import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk_
    import ttk
    import tkFileDialog
else:
    import tkinter as Tk_
    import tkinter.ttk as ttk
    import tkinter.filedialog as tkFileDialog


class PolyhedronViewer:
    """
    Window for viewing a hyperbolic polyhedron, either in the Poincare
    or Klein model.
    """

    def __init__(self, facedicts, root=None, title='Polyhedron Viewer',
                 container=None, bgcolor='#f4f4f4'):
        self.bgcolor = bgcolor
        self.style = SnapPyStyle()
        self.font = self.style.ttk_style.lookup('TLable', 'font')
        self.empty = (len(facedicts) == 0)
        self.title=title
        if root is None:
            if Tk_._default_root is None:
                root = Tk_.Tk()
            else:
                root = Tk_._default_root
            root.withdraw()
        self.root = root
        if container:
            self.window = window = container
        else:
            self.window = window = Tk_.Toplevel(master=root, class_='snappy')
            window.withdraw()
            window.title(title)
            window.protocol("WM_DELETE_WINDOW", self.close)
        self.menubar = None
        self.topframe = topframe = ttk.Frame(window)
        self.bottomframe = bottomframe = ttk.Frame(window)
        self.model_var=Tk_.StringVar(self.window, value='Klein')
        self.sphere_var=Tk_.IntVar(self.window, value=1)
        self.klein = ttk.Radiobutton(topframe, text='Klein',
                                     variable = self.model_var,
                                     value='Klein',
                                     command=self.new_model)
        self.poincare = ttk.Radiobutton(topframe, text='Poincaré',
                                        variable = self.model_var,
                                        value='Poincare',
                                        command=self.new_model)
        self.sphere = ttk.Checkbutton(topframe, text='',
                                      variable = self.sphere_var,
                                      command=self.new_model)
        self.spherelabel = spherelabel = Tk_.Text(topframe, height=1, width=3,
                                        relief=Tk_.FLAT, font=self.font,
                                        borderwidth=0, highlightthickness=0,
                                        background=bgcolor)
        spherelabel.tag_config("sub", offset=-4)
        spherelabel.insert(Tk_.END, 'S')
        spherelabel.insert(Tk_.END, '∞', 'sub')
        spherelabel.config(state=Tk_.DISABLED)
        if sys.platform == 'darwin':
            spherelabel.configure(background=self.style.groupBG)
        self.klein.grid(row=0, column=0, sticky=Tk_.W, padx=20, pady=(2,6))
        self.poincare.grid(row=0, column=1, sticky=Tk_.W, padx=20, pady=(2,6))
        self.sphere.grid(row=0, column=2, sticky=Tk_.W, padx=0, pady=(2,6))
        spherelabel.grid(row=0, column=3, sticky=Tk_.NW)
        topframe.pack(side=Tk_.TOP, fill=Tk_.X)
        self.widget = widget = OpenGLWidget(master=bottomframe,
                                            width=600,
                                            height=500,
                                            double=1,
                                            depth=1,
                                            help="""
Use mouse button 1 to rotate the polyhedron.

Releasing the button while moving will "throw" the polyhedron and make it keep spinning.

The slider controls zooming.  You will see inside the polyhedron if you zoom far enough.
""")
        widget.set_eyepoint(5.0)
        self.GL = GL_context()
        self.polyhedron = HyperbolicPolyhedron(facedicts,
                                               self.model_var,
                                               self.sphere_var)
        widget.redraw = self.polyhedron.draw
        widget.autospin_allowed = 1
        widget.set_background(.2, .2, .2)
        widget.grid(row=0, column=0, sticky=Tk_.NSEW)
        zoomframe = ttk.Frame(bottomframe)
        self.zoom = zoom = ttk.Scale(zoomframe, from_=100, to=0,
            length=500, orient=Tk_.VERTICAL, command=self.set_zoom)
        zoom.set(50)
        zoom.pack(side=Tk_.TOP, expand=Tk_.YES, fill=Tk_.Y)
        bottomframe.columnconfigure(0, weight=1)
        bottomframe.rowconfigure(0, weight=1)
        zoomframe.grid(row=0, column=1, sticky=Tk_.NS)
        bottomframe.pack(side=Tk_.TOP, expand=Tk_.YES, fill=Tk_.BOTH)
        self.build_menus()
        if container is None:
            if self.menubar:
                self.window.config(menu=self.menubar)
            window.deiconify()
        self.add_help()

    # Subclasses may override this, e.g. if there is a help menu already.
    def add_help(self):
        help = Tk_.Button(self.topframe, text = 'Help', width = 4,
                          borderwidth=0, highlightthickness=0,
                          background=self.bgcolor, command = self.widget.help)
        help.grid(row=0, column=4, sticky=Tk_.E, pady=3)
        self.topframe.columnconfigure(3, weight = 1)

    def export_stl(self):
        model = self.model_var.get()
        filename = tkFileDialog.asksaveasfilename(
            parent=self.window,
            title='Save %s model as STL file' % model,
            defaultextension = '.stl',
            filetypes = [
                ('STL files', '*.stl'),
                ('All files', '')])
        if filename == '':  # If user clicked cancel:
            return
        with open(filename, 'w') as output_file:
            n = 0
            for line in stl(self.polyhedron.facedicts, model=model.lower()):
                output_file.write(line)
                # This can take a long time so make sure the GUI stays alive.
                if n > 100:
                    self.root.update_idletasks()
                    n = 0

    def export_cutout_stl(self):
        model = self.model_var.get()
        filename = tkFileDialog.asksaveasfilename(
            parent=self.window,
            title='Save %s model cutout as STL file' % model,
            defaultextension = '.stl',
            filetypes = [
                ('STL files', '*.stl'),
                ('All files', '')])
        if filename == '':  # If user clicked cancel:
            return
        with open(filename, 'w') as output_file:
            n = 100
            for line in stl(self.polyhedron.facedicts, model=model.lower(), cutout=True):
                output_file.write(line)
                # This can take a long time so make sure the GUI stays alive.
                if n > 100:
                    self.root.update_idletasks()
                    n = 0

  # Subclasses may override this to provide menus.
    def build_menus(self):
        pass

  # Subclasses may override this to update menus, e.g. when embedded in a larger window.
    def update_menus(self, menubar):
        pass

    def close(self):
        self.polyhedron.destroy()
        self.window.destroy()

    def reopen(self):
        self.widget.tkRedraw()

    def reset(self):
        self.widget.autospin = 0
        self.widget.set_eyepoint(5.0)
        self.zoom.set(50)
        self.widget.tkRedraw()

    def set_zoom(self, x):
        t = float(x)/100.0
        self.widget.distance = t*1.0 + (1-t)*8.0
        self.widget.tkRedraw()

    def new_model(self):
        self.widget.tkRedraw()

    def new_polyhedron(self, new_facedicts):
        self.empty = (len(new_facedicts) == 0)
        self.polyhedron = HyperbolicPolyhedron(new_facedicts,
                                               self.model_var,
                                               self.sphere_var)
        self.widget.redraw = self.polyhedron.draw
        for n in range(5):
            self.widget.after(n*500, self.widget.tkRedraw)


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
