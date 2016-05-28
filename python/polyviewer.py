# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from snappy.CyOpenGL import *
import math
from time import time
try:
    import Tkinter as Tk_
    import ttk
    import tkFileDialog
except ImportError: #Python 3
    import tkinter as Tk_
    import tkinter.ttk


def facet_stl(vertex1, vertex2, vertex3):
    a = (vertex3[0]-vertex1[0], vertex3[1]-vertex1[1], vertex3[2]-vertex1[2])
    b = (vertex2[0]-vertex1[0], vertex2[1]-vertex1[1], vertex2[2]-vertex1[2])
    normal = (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])
    return ''.join([
        '  facet normal %f %f %f\n' % normal,
        '    outer loop\n',
        '      vertex %f %f %f\n' % tuple(vertex1),
        '      vertex %f %f %f\n' % tuple(vertex2),
        '      vertex %f %f %f\n' % tuple(vertex3),
        '    endloop\n',
        '  endfacet\n'
        ])

def tri_div(triangles):
    return [triangle for (x, y, z) in triangles for triangle in [
        (x, midpoint(x, y), midpoint(x, z)),
        (midpoint(y, x), y, midpoint(y, z)),
        (midpoint(z, x), midpoint(z, y), z),
        (midpoint(x, y), midpoint(y, z), midpoint(z, x))
        ]]

def midpoint(vertex1, vertex2):
    x1, y1, z1 = vertex1
    x2, y2, z2 = vertex2
    return ((x1+x2) / 2, (y1+y2) / 2, (z1+z2) / 2)

def projection(vertex, cutoff_radius=0.9):
    x, y, z = vertex
    D = x**2 + y**2 + z**2
    scale = 1 / (1 + math.sqrt(max(0, 1-D)))
    if scale >= cutoff_radius: scale = cutoff_radius
    return (scale*x, scale*y, scale*z)


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
        path = tkFileDialog.asksaveasfilename(
            parent=self.window,
            title='Save %s Model as STL file' % model,
            defaultextension = '.stl',
            filetypes = [
                ('STL files', '*.stl'),
                ('All files', '')])
        if path is None:
            return
        if model == 'Klein':
            output = self.klein_to_stl()
        elif model == 'Poincare':
            output = self.poincare_to_stl()
        else:
            raise ValueError('Unknown model')
        with open(path, 'w') as output_file:
            output_file.writelines(output)

    def export_cutout_stl(self):
        model = self.model_var.get()
        path = tkFileDialog.asksaveasfilename(
            parent=self.window,
            title='Save %s Model Cutout as STL file' % model,
            defaultextension = '.stl',
            filetypes = [
                ('STL files', '*.stl'),
                ('All files', '')])
        if path is None:
            return
        if model == 'Klein':
            output = self.klein_cutout()
        elif model == 'Poincare':
            output = self.poincare_cutout()
        else:
            raise ValueError('Unknown model')
        with open(path, 'w') as output_file:
            output_file.writelines(output)

    def klein_to_stl(self):
        output = ['solid\n']
        klein_faces = self.polyhedron.get_facedicts()
        for face in klein_faces:
            vertices = face['vertices']
            for i in range(len(vertices)-2):
                output.append(facet_stl(vertices[0], vertices[i+1], vertices[i+2]))
        output.append('endsolid')
        return output

    def poincare_to_stl(self, num_subdivisions=5):
        output = ['solid\n']
        klein_faces = self.polyhedron.get_facedicts()
        for face in klein_faces:
            vertices = face['vertices']
            for i in range(len(vertices)-2):
                triangles = [[vertices[0], vertices[i+1], vertices[i+2]]]
                for i in range(num_subdivisions):  # Subdivide.
                    triangles = tri_div(triangles)
                for triangle in triangles:
                    output.append(facet_stl(projection(triangle[0]), projection(triangle[1]), projection(triangle[2])))
        output.append('endsolid')
        return output

    def klein_cutout(self):
        output = ['solid\n']
        klein_faces = self.polyhedron.get_facedicts()
        for face in klein_faces:
            vertices = face['vertices']
            center = [sum(vertex[i] for vertex in vertices) / len(vertices) for i in range(3)]
            new_vertices = [[vertex[i] + (center[i] - vertex[i]) / 3 for i in range(3)] for vertex in vertices]
            new_inside_points = [[point[i] * 0.8 for i in range(3)] for point in new_vertices]
            for i in range(len(new_vertices)):
                output.append(facet_stl(new_vertices[i], new_inside_points[(i+1) % len(new_vertices)], new_inside_points[i]))
            for i in range(len(new_vertices)):
                output.append(facet_stl(new_vertices[i], new_vertices[(i+1) % len(new_vertices)], new_inside_points[(i+1) % len(new_vertices)]))
            for i in range(len(vertices)):
                output.append(facet_stl(vertices[i], new_vertices[(i+1) % len(vertices)], new_vertices[i]))
                output.append(facet_stl(vertices[i], vertices[(i+1) % len(vertices)], new_vertices[(i+1) % len(vertices)]))
                # We have to go in the opposite direction this time as the normal should point in towards O.
                point_list = [
                    (vertices[i], new_vertices[i], new_vertices[(i+1) % len(vertices)]),
                    (vertices[i], new_vertices[(i+1) % len(vertices)], vertices[(i+1) % len(vertices)])
                    ]
                scaled_point_list = [[[point[i] * 0.8 for i in range(3)] for point in triangle] for triangle in point_list]
                for triangle in scaled_point_list:
                    output.append(facet_stl(*triangle))
        output.append('endsolid')
        return output

    def poincare_cutout(self, num_subdivisions=3):
        start_time = time()
        output = ['solid\n']
        klein_faces = self.polyhedron.get_facedicts()
        for face in klein_faces:
            vertices = face['vertices']
            center = [sum(vertex[i] for vertex in vertices) / len(vertices) for i in range(3)]
            new_points = new_vertices = [[vertex[i] + (center[i] - vertex[i]) / 3 for i in range(3)] for vertex in vertices]
            # Subdivide new_points by inserting midpoints.
            for j in range(num_subdivisions):
                new_points = [point for point_midpoint in zip(new_points, [midpoint(new_points[i], new_points[(i+1) % len(new_points)]) for i in range(len(new_points))]) for point in point_midpoint]
            # Project and rescale to get the inside points.
            new_points = [projection(point) for point in new_points]
            new_inside_points = [[point[i] * 0.8 for i in range(3)] for point in new_points]
            for i in range(len(new_points)):
                output.append(facet_stl(new_points[i], new_inside_points[(i+1) % len(new_points)], new_inside_points[i]))
                output.append(facet_stl(new_points[i], new_points[(i+1) % len(new_points)], new_inside_points[(i+1) % len(new_points)]))
            for i in range(len(vertices)):
                triangles = [
                    (vertices[i], new_vertices[(i+1) % len(vertices)], new_vertices[i]),
                    (vertices[i], vertices[(i+1) % len(vertices)], new_vertices[(i+1) % len(vertices)])
                    ]
                for i in range(num_subdivisions):
                    triangles = tri_div(triangles)
                for triangle in triangles:
                    output.append(facet_stl(*[projection(triangle[i]) for i in range(3)]))
                # We have to go in the opposite direction this time as the normal should point in towards O.
                point_list = [[projection(vertex) for vertex in reversed(triangle)] for triangle in triangles]
                scaled_point_list = [[[point[i] * 0.8 for i in range(3)] for point in triangle] for triangle in point_list]
                for triangle in scaled_point_list:
                    output.append(facet_stl(*triangle))
        output.append('endsolid')
        print(time() - start_time)
        return output

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

