#! /usr/bin/env python
import os, sys, shutil, platform

this_python = sys.executable
this_pyinstaller = os.path.abspath(i
    os.path.join(this_python, '..', 'Scripts', 'pyinstaller'))

# The Inno Installer config file (*.iss) assumes a 64 bit binary.

if platform.architecture()[0] != '64bit' and '--32-bit' not in sys.argv:
    print("ERROR: Need to use a 64bit Python to build the apps")
    sys.exit(1)

try:
    import pyx
except ImportError:
    print("ERROR: Need to install PyX!")
    sys.exit(1)

try:
    import snappy_15_knots
except ImportError:
    print("ERROR: Need to install snappy_15_knots!")
    sys.exit(1)

if '--no-freshen' not in sys.argv:
    os.chdir("../windows_exe/../")
    os.system("git pull")
    os.system("rm dist/*.whl")
    os.system(this_python + " setup.py pip_install")
    os.chdir("windows_exe")

# Now build the .exe
os.system("rm -rf build dist InstallSnapPy.exe InstallSnapPy-Dbg.exe")
os.system(this_pyinstaller + " SnapPy_py3.spec")
os.system("iscc InnoSnapPy_py3.iss")
os.system(this_pyinstaller + " SnapPy_dbg.spec")
os.system("iscc InnoSnapPy_dbg.iss")
