#! /usr/bin/env python

import os, sys, re, glob

# We build the thing and install it twice to make sure the
# documentation is up to date.

python26 = "/Library/Frameworks/Python.framework/Versions/2.6/bin/python"
framework = '/Library/Frameworks/Python-10.5-intel.framework'
if not os.path.exists(framework):
    framework = '/Library/Frameworks/Python.framework'
print 'Using python from %s'%framework
python27 = os.path.join(framework, 'Versions', '2.7', 'bin', 'python')

os.chdir("../")
os.system("hg pull")
os.system("hg up")
os.system(python27 + " setup.py clean")
for python in [python26, python27]:
    os.system(python + " setup.py install")
    os.system(python + " setup.py build_docs install")
    
# Now build the .app

os.chdir("SnapPyApp")
os.system(python27 + " setup.py py2app")

# Make things a little smaller.

os.system("rm -rf dist/SnapPy.app/Contents/Frameworks/Tcl.framework/Versions/*/Resources/English.lproj/ActiveTcl-*")
os.system("rm -rf dist/SnapPy.app/Contents/Frameworks/Tk.framework/Versions/*/Resources/Scripts/demos")

# Make sure we use the correct version of Tk (8.5 aat the moment):
os.system("pushd dist/SnapPy.app/Contents/Frameworks/Tcl.framework/Versions/; "
          "rm Current; ln -s 8.5 Current; popd")
os.system("pushd dist/SnapPy.app/Contents/Frameworks/Tk.framework/Versions/; "
          "rm Current; ln -s 8.5 Current; popd")

# Then make the disk image file.  

os.chdir("dmg-maker")
os.system("./dmg-maker.py")

# Now put it on the webpage:

user = os.environ['USER']
if user in ['nmd', 'dunfield']:
    print "Hi there Nathan..."
    address = "nmd@shell.math.uic.edu"
if user == 'culler':
    print "Hi there Marc..."
    address = "culler@threlfall.math.uic.edu"


for file in glob.glob("../../dist/*-intel.egg"):
    copy = file.replace("-intel", "-fat")
    os.system("cp " + file + " " + copy)
    
os.system("chmod g+w SnapPy.dmg ../../dist/*.egg")
os.system("scp -p SnapPy.dmg ../../dist/*.egg %s:/afs/math.uic.edu/www/t3m/SnapPy-nest" % address)
