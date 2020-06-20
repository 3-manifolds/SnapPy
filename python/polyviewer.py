# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import time

from .gui import *
from .CyOpenGL import *
from .export_stl import stl
from . import filedialog
from plink.ipython_tools import IPythonTkRoot

class PolyhedronViewer(ttk.Frame):
    """
    Displays a hyperbolic polyhedron, either in the Poincare or Klein model.
    """

    def __init__(self, master, facedicts=[], **kwargs):
        ttk.Frame.__init__(self, master)
        self.empty = (len(facedicts) == 0)
        self.style = style = SnapPyStyle()
        self.bgcolor = kwargs.get('bgcolor', self.style.windowBG)
        self.menubar = None
        self.main_window = kwargs.get('main_window', None)
        self.topframe = topframe = ttk.Frame(self)
        self.bottomframe = bottomframe = ttk.Frame(self)
        self.model_var=Tk_.StringVar(self, value='Klein')
        self.sphere_var=Tk_.IntVar(self, value=1)
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
                                                  relief=Tk_.FLAT,
                                                  font=self.style.font,
                                                  borderwidth=0,
                                                  highlightthickness=0,
                                                  background=self.bgcolor)
        spherelabel.tag_config("sub", offset=-4)
        spherelabel.insert(Tk_.END, 'S')
        spherelabel.insert(Tk_.END, '∞', 'sub')
        spherelabel.config(state=Tk_.DISABLED)
        if sys.platform == 'darwin':
            if master:
                spherelabel.configure(background=self.style.groupBG)
            else:
                spherelabel.configure(background=self.style.windowBG)
        self.klein.grid(row=0, column=0, sticky=Tk_.W, padx=20, pady=(2,6))
        self.poincare.grid(row=0, column=1, sticky=Tk_.W, padx=20, pady=(2,6))
        self.sphere.grid(row=0, column=2, sticky=Tk_.W, padx=0, pady=(2,6))
        spherelabel.grid(row=0, column=3, sticky=Tk_.NW)
        topframe.pack(side=Tk_.TOP, fill=Tk_.X)
        self.widget = widget = OpenGLPerspectiveWidget(master = bottomframe,
                                                       width = 600,
                                                       height = 500,
                                                       double = 1,
                                                       depth = 1,
                                                       help = """
Use mouse button 1 to rotate the polyhedron.

Releasing the button while moving will "throw" the polyhedron and make it keep spinning.

The slider controls zooming.  You will see inside the polyhedron if you zoom far enough.
""")
        widget.set_eyepoint(5.0)
        cyglSetStandardLighting()
        self.polyhedron = HyperbolicPolyhedron(facedicts, self.model_var,
            self.sphere_var, togl_widget=self.widget)
        widget.redraw_impl = self.polyhedron.draw
        widget.autospin_allowed = 1
        widget.set_background(.2, .2, .2)
        widget.grid(row=0, column=0, sticky=Tk_.NSEW)
        zoomframe = ttk.Frame(bottomframe)
        self.zoom = zoom = ttk.Scale(zoomframe, from_=0, to=100,
            orient=Tk_.VERTICAL, command=self.set_zoom)
        zoom.set(50)
        zoom.pack(side=Tk_.TOP, expand=Tk_.YES, fill=Tk_.Y)
        bottomframe.columnconfigure(0, weight=1)
        bottomframe.rowconfigure(0, weight=1)
        zoomframe.grid(row=0, column=1, sticky=Tk_.NS)
        bottomframe.pack(side=Tk_.TOP, expand=Tk_.YES, fill=Tk_.BOTH)
        self.build_menus()
        if isinstance(master, Tk_.Toplevel) and self.menubar:
            master.config(menu=self.menubar)
        self.add_help()
        # Added to avoid occasional missing faces in the browser.
        self.update_idletasks()

    def apply_prefs(self, prefs):
        # We don't have any preferences to apply
        return

    # Subclasses may override this, e.g. if there is a help menu already.
    def add_help(self):
        help = ttk.Button(self.topframe, text = 'Help', width = 4,
                          command = self.widget.help)
        help.grid(row=0, column=4, sticky=Tk_.E, padx=18)
        self.topframe.columnconfigure(3, weight = 1)

    def export_stl(self):
        model = self.model_var.get()
        file = filedialog.asksaveasfile(
            parent=self.master,
            title='Save %s model as STL file' % model,
            defaultextension = '.stl',
            filetypes = [
                ('STL files', '*.stl'),
                ('All files', '')])
        if file:
            n = 0
            for line in stl(self.polyhedron.facedicts, model=model.lower()):
                file.write(line)
                # This can take a long time so make sure the GUI stays alive.
                if n > 100:
                    self.root.update_idletasks()
                    n = 0
            file.close()

    def export_cutout_stl(self):
        model = self.model_var.get()
        file = filedialog.asksaveasfile(
            parent=self.master,
            title='Save %s model cutout as STL file' % model,
            defaultextension = '.stl',
            filetypes = [
                ('STL files', '*.stl'),
                ('All files', '')])
        if file:
            n = 100
            for line in stl(self.polyhedron.facedicts, model=model.lower(), cutout=True):
                file.write(line)
                # This can take a long time so make sure the GUI stays alive.
                if n > 100:
                    self.root.update_idletasks()
                    n = 0
            file.close()

  # Subclasses may override this to provide menus.
    def build_menus(self):
        pass

  # Subclasses may override this to update menus, e.g. when embedded in a larger window.
    def update_menus(self, menubar):
        pass

    def redraw(self):
        self.widget.redraw_if_initialized()

    def reset(self):
        self.widget.autospin = 0
        self.widget.set_eyepoint(5.0)
        self.zoom.set(50)
        self.widget.redraw_if_initialized()

    def set_zoom(self, x):
        t = (100.0-float(x))/100.0
        self.widget.distance = t*1.0 + (1-t)*8.0
        self.widget.redraw_if_initialized()

    def new_model(self):
        self.widget.redraw_if_initialized()

    def new_polyhedron(self, new_facedicts):
        self.empty = (len(new_facedicts) == 0)
        self.widget.tk.call(self.widget._w, 'makecurrent')
        try:
            self.polyhedron.delete_resource()
        except AttributeError:
            pass
        self.polyhedron = HyperbolicPolyhedron(new_facedicts, self.model_var,
            self.sphere_var, togl_widget=self.widget)
        self.widget.redraw_impl = self.polyhedron.draw
        self.widget.redraw_if_initialized()

    def delete_resource(self):
        try:
            self.polyhedron.delete_resource()
        except AttributeError:
            pass

    def test(self):
        X = 100
        self.widget.event_generate('<Button-1>', x=X, y=300, warp=True)
        self.update_idletasks()
        for n in range(10):
            X += 30
            time.sleep(0.1)
            self.widget.event_generate('<B1-Motion>', x=X, y=300, warp=True)
        self.widget.event_generate('<ButtonRelease-1>', x=X+30, y=300, warp=True)
        self.update_idletasks()
        time.sleep(0.5)
        self.model_var.set('Poincare')
        self.update_idletasks()
        self.new_model()
        self.update_idletasks()
        time.sleep(1.0)
        self.model_var.set('Klein')
        self.update_idletasks()
        self.new_model()
        self.update_idletasks()
        time.sleep(0.5)

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
    PV = ViewerWindow(PolyhedronViewer, testpoly)
    PV.mainloop()
