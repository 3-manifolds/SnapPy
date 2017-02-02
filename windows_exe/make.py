#! /usr/bin/env python
from __future__ import print_function

import os, sys, re, shutil

this_python = sys.executable
this_pyinstaller = os.path.abspath(
    os.path.join(this_python, '..', 'Scripts', 'pyinstaller'))

try:
    import pyx
except ImportError:
    print("ERROR: Need to install PyX!")
    sys.exit()

os.chdir("../windows_exe/../")
os.system("hg pull")
os.system("hg update")
os.system("rm dist/*.egg")

os.system(this_python + " setup.py install")
os.system(this_python + " setup.py build_docs")
os.system(this_python + " setup.py install")

# Now build the .exe

os.chdir("windows_exe")
os.system("rm -rf build dist InstallSnappy.exe")
if sys.version_info.major == 2:
    os.system(this_pyinstaller + " SnapPy_py2.spec")
else:
    os.system(this_pyinstaller + " SnapPy_py3.spec")
    
print("Starting the app to force lib2to3 to build pickles.")
print("Close the app to continue.")
os.system(os.path.join("dist", "SnapPy", "SnapPy.exe"))

# Build the Inno Setup installer
os.system("iscc InnoSnapPy.iss")
