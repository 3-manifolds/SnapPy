#! /usr/bin/env python

import os, sys, re

# We build the thing and install it twice to make sure the
# documentation is up to date.

os.chdir("../")
os.system("python setup.py install")
os.chdir("doc-source")
os.system("make install")
os.chdir("../")
os.system("python setup.py install")

# Now build the .app

os.chdir("SnapPyApp")
os.system("python setup.py clean py2app")

# Make things a little smaller.

os.system("rm -rf dist/SnapPy.app/Contents/Frameworks/Tcl.framework/Versions/8.4/Resources/English.lproj/ActiveTcl-8.4")
os.system("rm -rf dist/SnapPy.app/Contents/Frameworks/Tk.framework/Versions/8.4/Resources/Scripts/demos")

# Then the disk image file.  

os.chdir("dmg-maker")
os.system("./dmg-maker.py")

# Now put it on the webpage:

user = os.environ['USER']
if user in ['nmd', 'dunfield']:
    print "Hi there Nathan..."
    os.system("scp SnapPy.dmg t3m@shell.math.uic.edu:public_html/")
if user == 'culler':
    print "Hi there Marc..."
    os.system("scp SnapPy.dmg culler@shell.math.uic.edu:~t3m/public_html/")
