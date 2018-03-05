#! /usr/bin/env python
from __future__ import print_function

import os, sys, shutil, platform

this_python = sys.executable
this_pyinstaller = os.path.abspath(
    os.path.join(this_python, '..', 'Scripts', 'pyinstaller'))

# We currently build the Windows apps as 32 bit apps only.  The reason for
# this is a mysterious unresolved issue with the 64 bit app, namely
# that when installed in C:\Program Files it takes close to one minute
# to start up.  While it does start up normally if installed elsewhere,
# we want to resolve this issue before releasing a 64 bit app.

if platform.architecture()[0] != '32bit':
    print("ERROR: Need to use a 32bit Python to build the apps")
    sys.exit()

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
os.system("rm -rf build dist InstallSnappy-Python%s.exe"%sys.version_info.major)
if sys.version_info.major == 2:
    os.system(this_pyinstaller + " SnapPy_py2.spec")
    os.system("iscc InnoSnapPy_py2.iss")
else:
    os.system(this_pyinstaller + " SnapPy_py3.spec")
    os.system("iscc InnoSnapPy_py3.iss")
    os.system(this_pyinstaller + " SnapPy_dbg.spec")
    os.system("iscc InnoSnapPy_dbg.iss")
