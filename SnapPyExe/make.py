#! /usr/bin/env python

import os, sys, re

# We build the thing and install it twice to make sure the
# documentation is up to date.
glut_dll = '/c/python26/Lib/site-packages/pyopengl-3.0.0-py2.6.egg/OpenGL/DLLS/glut32.dll'


os.chdir("../")
os.system("python setup.py install")
os.chdir("doc-source")
os.system("make install")
os.chdir("../")
os.system("python setup.py install")

# Now build the .app

os.chdir("SnapPyExe")
os.system("rm -rf build InstallSnappy.exe SnapPy")
os.system('cp %s .'%glut_dll)
os.system("python setup.py py2exe")

# Make things a little smaller.

os.system("rm -rf SnapPy/tcl/tk8.5/demos")

# Build the Inno Setup installer
os.system("compil32 /cc InnoSnapPy.iss")
