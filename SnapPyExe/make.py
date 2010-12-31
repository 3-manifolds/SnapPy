#! /usr/bin/env python

import os, sys, re

python26 = "c:\Python26\python.exe "
python27 = "c:\Python27\python.exe " 

os.chdir("../SnapPyExe/../")
os.system("hg pull")
os.system("hg update")

#for python in [python27]:
for python in [python26, python27]:
    os.system(python + "setup.py build -c mingw32")
    os.system(python + "setup.py install")
    os.system(python + "setup.py build_docs")
    os.system(python + "setup.py install")

# Now build the .exe

os.chdir("SnapPyExe")
os.system("rm -rf build InstallSnappy.exe SnapPy")
os.system(python27 + "setup.py py2exe")

# Make things a little smaller.

os.system("rm -rf SnapPy/tcl/tk8.5/demos")

# Build the Inno Setup installer

os.system("compil32 /cc InnoSnapPy.iss")

# Copy the installer to the website

address = "t3m@shell.math.uic.edu"
os.system("scp InstallSnapPy.exe ../dist/*.egg %s:/home/www/t3m/public_html/SnapPy-nest" % address)
