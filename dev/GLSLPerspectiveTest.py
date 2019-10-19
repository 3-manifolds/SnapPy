from __future__ import print_function

import tkinter as Tk_
import tkinter.ttk as ttk

from snappy.CyOpenGL import *
from math import sqrt

vertex_shader_source = b"""
#version 150

// Note that GLSL ES 3.0 is based of GLSL 3.30 and is used for WebGL 2.0.
// GLSL 1.50 came with OpenGL 3.2.

in vec4 position;
uniform mat4 MVPMatrix;
uniform mat4 ModelViewMatrix;
uniform mat4 ProjectionMatrix;

void main()
{
    gl_Position = MVPMatrix*position;
}

"""

fragment_shader_source = b"""
#version 150

//varying float myRed;

out vec4 out_FragColor;

void main()
{
  if (gl_FrontFacing)
    {
      out_FragColor = vec4(1.0,0.0,0.0,1.0);;
    }
  else
    {
      out_FragColor = vec4(0.0,0.0,1.0,1.0);
    }}
"""

class PerspectiveTester:
    def __init__(self, window=None):
        if window == None:
            window = Tk_.toplevel(Tk_._default_root)
        self.window = window
        window.grid_rowconfigure(0, weight=1)
        window.grid_columnconfigure(0, weight=1)
        window.title('GLSL Perspective Test')
        self.widget = widget = GLSLPerspectiveWidget(window,
            width=600, height=500, double=1, depth=1)
        widget.grid(row=0, column=0, sticky=Tk_.NSEW)
        zoomer = ttk.Scale(window, from_=100, to=0,
                           orient=Tk_.VERTICAL,
                           command=self.set_zoom)
        zoomer.set(50)
        zoomer.grid(row=0, column=1, sticky=Tk_.NSEW)
        window.bind("<Left>", lambda event: self.rotate(-0.1, 0, 1.0, 0))
        window.bind("<Right>", lambda event: self.rotate(0.1, 0, 1.0, 0))
        window.bind("<Up>", lambda event: self.rotate(-0.1, 1.0, 0, 0))
        window.bind("<Down>", lambda event: self.rotate(0.1, 1.0, 0, 0))
        widget.bind("<Button-1>", self._mouse_press)
        widget.bind("<B1-Motion>", self._mouse_drag)
        widget.make_current()
        # Objects cannot be created without a GL context.  We can do it now.
        widget.add_object(
            Triangle(
                widget,
                GLSLProgram(vertex_shader_source, fragment_shader_source),
                (( -1.0, -0.867, 0.0),
                 (  1.0, -0.867, 0.0),
                 ( -0.0,  0.867, 0.0))))
        widget.redraw_if_initialized()
        # Print out the GL version for testing.
        print(get_gl_string('GL_VERSION'))

    def set_zoom(self, x):
        t = float(x)/100.0
        self.widget.distance = t*1.0 + (1-t)*8.0
        self.widget.redraw_if_initialized()

    def rotate(self, angle, x, y, z):
        self.widget.rotate(angle, x, y, z)
        self.widget.redraw_if_initialized()

    def _mouse_press(self, event):
        """
        Record the mouse position when button 1 is pressed.
        """
        self.mouse_x = event.x
        self.mouse_y = event.y

    def _mouse_drag(self, event):
        """
        Rotate the model around an axis perpendicular to the mouse displacement by
        an angle proportional to the distance the mouse was moved.
        """
        dx = float(event.x - self.mouse_x)
        dy = float(event.y - self.mouse_y)
        self.mouse_x = event.x
        self.mouse_y = event.y
        norm = sqrt(dx*dx + dy*dy)
        if norm > .000001:
            self.rotate(-0.01*norm, -dy/norm, dx/norm, 0)

def main():
    root = Tk_.Tk()
    tester = PerspectiveTester(root)
    root.mainloop()

if __name__ == '__main__':
    main()
