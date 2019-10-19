from __future__ import print_function

import tkinter as Tk_
from snappy.CyOpenGL import *

fragment_shader_source = b"""
#version 150

out vec4 out_FragColor;
uniform float ViewportWidth;
uniform float ViewportHeight;

void main()
{
    float x = gl_FragCoord[0], y = gl_FragCoord[1];
    out_FragColor = vec4(1.0, x/ViewportWidth, y/ViewportHeight, 1.0);
}
"""  

class ImageShaderWidget(SimpleImageShaderWidget):
    def get_uniform_bindings(self, width, height):
        return {
            'ViewportWidth': ('float', width),
            'ViewportHeight': ('float', height) }

def create_widget(toplevel):
    widget = ImageShaderWidget(toplevel, fragment_shader_source,
        width=600, height=500, double=1, depth=1)
    widget.make_current()
    print(get_gl_string('GL_VERSION'))
    toplevel.grid_rowconfigure(0, weight=1)
    toplevel.grid_columnconfigure(0, weight=1)
    widget.grid(row = 0, column = 0, sticky = Tk_.NSEW)
    return widget

def main():
    root = Tk_.Tk()
    root.title('Image Shader Test')
    widget = create_widget(root)
    root.mainloop()
    
if __name__ == '__main__':
    main()
