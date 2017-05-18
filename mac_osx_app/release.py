#! /usr/bin/env python
import os, sys, re, glob

# We build the thing and install it twice to make sure the
# documentation is up to date.

framework = '/Library/Frameworks/Python.framework'
print ('Using python from %s'%framework)
if sys.version_info.major == 2:
    python = os.path.join(framework, 'Versions', '2.7', 'bin', 'python')
    tk_ver = os.popen(python + ' -c "import Tkinter; print(Tkinter.TkVersion)"').read().strip()
else:
    python = os.path.join(framework, 'Versions', '3.6', 'bin', 'python3')
    tk_ver = os.popen(python + ' -c "import tkinter; print(tkinter.TkVersion)"').read().strip()


os.chdir("../")
os.system("hg pull")
os.system("hg up")
os.system(python + " setup.py install")
os.system(python + " setup.py build_docs install")
    
# Now build the .app

os.chdir("mac_osx_app")
os.system(python + " setup.py py2app")

# We have to specify the Tcl and Tk frameworks for python 3.6, but then
# we get all versions of the framework which happen to be present on the
# build system.  We remove the ones we don't need.
tk_versions = "dist/SnapPy.app/Contents/Frameworks/Tk.framework/Versions"
tcl_versions = "dist/SnapPy.app/Contents/Frameworks/Tcl.framework/Versions"
for version in os.listdir(tk_versions):
    if version != tk_ver and version != 'Current':
        os.system("rm -rf %s"%os.path.join(tk_versions, version)) 
        os.system("rm -rf %s"%os.path.join(tcl_versions, version)) 

# Make things a little smaller.

os.system("rm -rf dist/SnapPy.app/Contents/Frameworks/Tcl.framework/Versions/*/Resources/English.lproj/ActiveTcl-*")
os.system("rm -rf dist/SnapPy.app/Contents/Frameworks/Tk.framework/Versions/*/Resources/Scripts/demos")

# Make sure we use the correct version of Tk:

os.system("pushd dist/SnapPy.app/Contents/Frameworks/Tcl.framework/Versions/; "
          "rm Current; ln -s " + tk_ver + " Current; popd")
os.system("pushd dist/SnapPy.app/Contents/Frameworks/Tk.framework/Versions/; "
          "rm Current; ln -s " + tk_ver + " Current; popd")

# Add a symlink so that py2app 0.13 can find the Tcl Scripts directory
os.system("pushd dist/SnapPy.app/Contents/; mkdir lib; cd lib; ln -s ../Frameworks/Tk.Framework/Versions/Current/Resources/Scripts tk%s; popd"%tk_ver)

# Then make the disk image file.  

os.chdir("dmg-maker")
os.system(python + " ./dmg-maker.py")
