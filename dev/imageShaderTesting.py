from __future__ import print_function

import tkinter as Tk_
from snappy.CyOpenGL import *

vertex_shader_source = b"""
#version 150

// Note that GLSL ES 3.0 is based of GLSL 3.30 and is used for WebGL 2.0.
// GLSL 1.50 came with OpenGL 3.2.

//varying float myRed;

in vec4 position;

void main()
{
    gl_Position = position; // gl_Vertex;//gl_ProjectionMatrix * gl_ModelViewMatrix * gl_Vertex;
  //  myRed = gl_Position.x;
}

"""

fragment_shader_source = b"""
#version 150

//varying float myRed;

/*
uniform MatrixBlock
{
  mat4 projection[];
  mat4 modelview[];
};
*/

out vec4 out_FragColor;

void main()
{
    out_FragColor = vec4(1.0,0.0,0.0,1.0);
}
"""  

def create_widget(toplevel):

    # imageShader = ImageShader()

    #widget = TestImageShaderOpenGLWidget(
    #    imageShader = imageShader,
    #    master = toplevel, width = 600, height = 500, double = 1, depth = 1)

    widget = SimpleImageShaderOpenGLWidget(
        vertex_shader_source,
        fragment_shader_source,
        master = toplevel,
        width = 600, height = 500, double = 1, depth = 1)
    widget.make_current()
    print(get_gl_string('GL_VERSION'))

    #imageShader.add_source(
#        b"""
#""",    

    widget.grid(row = 0, column = 0, sticky = Tk_.NSEW)
    return widget

def main():
    root = Tk_.Tk()
    widget = create_widget(root)
    root.mainloop()
    
if __name__ == '__main__':
    main()
