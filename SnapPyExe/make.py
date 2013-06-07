#! /usr/bin/env python

import os, sys, re

python26 = "c:\Python26\python.exe "
python27 = "c:\Python27\python.exe " 

os.chdir("../SnapPyExe/../")
os.system("hg pull")
os.system("hg update")
os.system("rm dist/*.egg")

for python in [python27]:
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

address = "nmd@shell.math.uic.edu"
raw_input('Hit any key when ready to begin copying to t3m:')
os.system("scp InstallSnapPy.exe ../dist/*.egg %s:/afs/math.uic.edu/www/t3m/SnapPy-nest" % address)
