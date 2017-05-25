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
os.system("rm -rf build dist InstallSnappy-Python2.exe InstallSnappy-Python3.exe")
if sys.version_info.major == 2:
    os.system(this_pyinstaller + " SnapPy_py2.spec")
    os.system("iscc InnoSnapPy_py2.iss")
else:
    os.system(this_pyinstaller + " SnapPy_py3.spec")
    os.system("iscc InnoSnapPy_py3.iss")
    os.system(this_pyinstaller + " SnapPy_dbg.spec")
    os.system("iscc InnoSnapPy_dbg.iss")
