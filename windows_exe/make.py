#! /usr/bin/env python

import os, sys, re

python27 = "c:\Python27\python.exe "

try:
    import pyx
except ImportError:
    print "ERROR: Need to install PyX!"
    sys.exit()

os.chdir("../windows_exe/../")
os.system("hg pull")
os.system("hg update")
os.system("rm dist/*.egg")

for python in [python27]:
    os.system(python + "setup.py install")
    os.system(python + "setup.py build_docs")
    os.system(python + "setup.py install")

# Now build the .exe

os.chdir("windows_exe")
os.system("rm -rf build InstallSnappy.exe SnapPy")
os.system(python27 + "setup.py py2exe")

# Make things a little smaller.

os.system("rm -rf SnapPy/tcl/tk8.5/demos")

# Build the Inno Setup installer

os.system("compil32 /cc InnoSnapPy.iss")

# Copy the installer to the website

address = "nmd@shell.math.uic.edu"
raw_input('Hit any key when ready to begin copying to t3m:')
os.system("scp InstallSnapPy.exe %s:/afs/math.uic.edu/www/t3m/SnapPy-nest" % address)
