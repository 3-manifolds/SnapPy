#! /usr/bin/env python

import os, sys, re

glut_dll = '/c/Python27/Lib/site-packages/pyopengl-3.0.1-py2.7-win32.egg/OpenGL/DLLS/glut32.dll'

# First install an egg with no docs

os.chdir("../SnapPyExe/../")
os.system("hg pull")
os.system("hg update")
os.system("python setup.py install")

# then go back and build the docs

os.system("python setup.py build_docs")
os.system("python setup.py install")

# Now build the .exe

os.chdir("SnapPyExe")
os.system("rm -rf build InstallSnappy.exe SnapPy")
os.system('cp %s .'%glut_dll)
os.system("python setup.py py2exe")

# Make things a little smaller.

os.system("rm -rf SnapPy/tcl/tk8.5/demos")

# Build the Inno Setup installer

os.system("compil32 /cc InnoSnapPy.iss")

# Copy the installer to the website

address = "t3m@shell.math.uic.edu"
os.system("scp InstallSnapPy.exe %s:/home/www/t3m/public_html/SnapPy-nest" % address)
os.system("scp ../../dist/*.egg %s:/home/www/t3m/public_html/SnapPy-nest" % address)
