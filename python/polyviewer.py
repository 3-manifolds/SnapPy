# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from snappy.CyOpenGL import *
from .export_stl import stl
try:
    import Tkinter as Tk_
    import ttk
    import tkFileDialog
except ImportError: #Python 3
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
        self.font = ttk.Style().lookup('TLable', 'font')
        self.empty = (len(facedicts) == 0)
        self.title=title
        if root is None:
            if Tk_._default_root is None:
                root = Tk_.Tk()
                root.iconify()
            else:
                root = Tk_._default_root
        self.root = root
        if container:
            self.window = window = container
        else:
            self.window = window = Tk_.Toplevel(master=root, class_='snappy')
            window.withdraw()
            window.title(title)
            window.protocol("WM_DELETE_WINDOW", self.close)
        self.menubar = None
        self.topframe = topframe = Tk_.Frame(window, borderwidth=0,
                                             relief=Tk_.FLAT, background=bgcolor)
        self.bottomframe = bottomframe = Tk_.Frame(window, borderwidth=0,
                                             relief=Tk_.FLAT)
        self.model_var=Tk_.StringVar(value='Klein')
        self.sphere_var=Tk_.IntVar(value=1)
        radiobutton_options = {
            'command' : self.new_model,
            'background' : bgcolor,
            'activebackground' : bgcolor,
            'highlightthickness' : 0,
            'borderwidth' : 0}
        self.klein = Tk_.Radiobutton(topframe, text='Klein',
                                     variable = self.model_var,
                                     value='Klein',
                                     **radiobutton_options)
        self.poincare = Tk_.Radiobutton(topframe, text='Poincaré',
                                        variable = self.model_var,
                                        value='Poincare',
                                        **radiobutton_options)
        self.sphere = Tk_.Checkbutton(topframe, text='',
                                      variable = self.sphere_var,
                                      **radiobutton_options)
        self.spherelabel = Tk_.Text(topframe, height=1, width=3,
                                    relief=Tk_.FLAT, font=self.font,
                                    borderwidth=0, highlightthickness=0,
                                    background=bgcolor)
        self.spherelabel.tag_config("sub", offset=-4)
        self.spherelabel.insert(Tk_.END, 'S')
        self.spherelabel.insert(Tk_.END, '∞', 'sub')
        self.spherelabel.config(state=Tk_.DISABLED)
        self.klein.grid(row=0, column=0, sticky=Tk_.W, padx=20, pady=(2,6))
        self.poincare.grid(row=0, column=1, sticky=Tk_.W, padx=20, pady=(2,6))
        self.sphere.grid(row=0, column=2, sticky=Tk_.W, padx=0, pady=(2,6))
        self.spherelabel.grid(row=0, column=3, sticky=Tk_.NW)
        topframe.pack(side=Tk_.TOP, fill=Tk_.X)
        self.widget = widget = OpenGLWidget(master=bottomframe,
                                            width=809,
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
        zoomframe = Tk_.Frame(bottomframe, borderwidth=0, relief=Tk_.FLAT,
                              background=self.bgcolor)
        self.zoom = zoom = Tk_.Scale(zoomframe, showvalue=0, from_=100, to=0,
                                     command=self.set_zoom, width=11,
                                     troughcolor=self.bgcolor, borderwidth=1,
                                     relief=Tk_.FLAT)
        zoom.set(50)
        spacer = Tk_.Frame(zoomframe, height=14, borderwidth=0, relief=Tk_.FLAT,
                           background=self.bgcolor)
        zoom.pack(side=Tk_.TOP, expand=Tk_.YES, fill=Tk_.Y)
        spacer.pack()
        bottomframe.columnconfigure(0, weight=1)
        bottomframe.rowconfigure(0, weight=1)
        zoomframe.grid(row=0, column=1, sticky=Tk_.NS)
        bottomframe.pack(side=Tk_.TOP, expand=Tk_.YES, fill=Tk_.BOTH)
        self.build_menus()
        if container is None:
            if self.menubar:
                self.window.config(menu=self.menubar)
            window.deiconify()
            window.update() # Seems to avoid a race condition with togl
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
            for line in stl(self.polyhedron.facedicts, model=model.lower()):
                self.root.update()  # This can take a long time so make sure the GUI stays alive.
                output_file.write(line)

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
            for line in stl(self.polyhedron.facedicts, model=model.lower(), cutout=True):
                self.root.update()  # This can take a long time so make sure the GUI stays alive.
                output_file.write(line)

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

