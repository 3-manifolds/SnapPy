#! /usr/bin/env python

import os, sys, re, glob

# We build the thing and install it twice to make sure the
# documentation is up to date.

framework = '/Library/Frameworks/Python.framework'
print 'Using python from %s'%framework
python27 = os.path.join(framework, 'Versions', '2.7', 'bin', 'python')

os.chdir("../")
os.system("hg pull")
os.system("hg up")
os.system(python27 + " setup.py install")
os.system(python27 + " setup.py build_docs install")
    
# Now build the .app

os.chdir("mac_osx_app")
os.system(python27 + " setup.py py2app")

# Make things a little smaller.

os.system("rm -rf dist/SnapPy.app/Contents/Frameworks/Tcl.framework/Versions/*/Resources/English.lproj/ActiveTcl-*")
os.system("rm -rf dist/SnapPy.app/Contents/Frameworks/Tk.framework/Versions/*/Resources/Scripts/demos")

tk_ver = os.popen(python27 + ' -c "import Tkinter; print(Tkinter.TkVersion)"').read().strip()

# Make sure we use the correct version of Tk:

os.system("pushd dist/SnapPy.app/Contents/Frameworks/Tcl.framework/Versions/; "
          "rm Current; ln -s " + tk_ver + " Current; popd")
os.system("pushd dist/SnapPy.app/Contents/Frameworks/Tk.framework/Versions/; "
          "rm Current; ln -s " + tk_ver + " Current; popd")

# Add a symlink so that py2app 0.10 can find the Tcl Scripts directory
os.system("pushd dist/SnapPy.app/Contents/; mkdir lib; cd lib; ln -s ../Frameworks/Tk.Framework tk"
          + tk_ver + "; popd")

# Then make the disk image file.  

os.chdir("dmg-maker")
os.system("./dmg-maker.py")
