#! /usr/bin/env python

import os, sys, re, shutil

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
os.system("rm -rf build dist InstallSnappy.exe")
os.system("pyinstaller SnapPy.spec")

print "Starting the app to force lib2to3 to build pickles."
print "Close the app to continue."
os.system(os.path.join("dist", "SnapPy", "SnapPy.exe"))

# Build the Inno Setup installer
os.system("iscc InnoSnapPy.iss")
