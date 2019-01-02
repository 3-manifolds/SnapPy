#! /usr/bin/env python
import os, sys, re
from glob import glob
from subprocess import call, Popen, PIPE
from math import ceil

def get_tk_ver(python):
    """
    Figure out which version of Tk is used by this python.
    """
    out, errors = Popen([python, "-c", "import _tkinter; print(_tkinter.TK_VERSION)"],
                    stdout = PIPE).communicate()
    return out.strip()

def freshen_SnapPy(python):
    """
    Build SnapPy and install it twice to make sure the documentation is
    up to date.
    """
    os.chdir("../")
    call(["hg", "pull"])
    call(["hg", "up"])
    call([python, "setup.py", "pip_install"])
    call([python, "setup.py", "build_docs", "pip_install"])
    os.chdir("mac_osx_app")
    
def build_app(python):
    """
    Build the standalone app bundle.
    """
    tk_ver = get_tk_ver(python)
    call([python, "setup.py", "py2app"])
    # We have to specify the Tcl and Tk frameworks for python 3.6, but then
    # we get all versions of the framework which happen to be present on the
    # build system.  We remove the ones we don't need.
    tk_versions = "dist/SnapPy.app/Contents/Frameworks/Tk.framework/Versions"
    tcl_versions = "dist/SnapPy.app/Contents/Frameworks/Tcl.framework/Versions"
    for version in os.listdir(tk_versions):
        if version != tk_ver and version != 'Current':
            os.system("rm -rf %s"%os.path.join(tk_versions, version)) 
            os.system("rm -rf %s"%os.path.join(tcl_versions, version)) 

def cleanup_app(python):
    """
    Tidy things up.
    """
    tk_ver = get_tk_ver(python)
    # Make things a little smaller.
    framework_dir = "dist/SnapPy.app/Contents/Frameworks/"
    junk = glob(framework_dir + "Tcl.framework/Versions/*/Resources/English.lproj/ActiveTcl-*")
    junk += glob(framework_dir + "Tk.framework/Versions/*/Resources/Scripts/demos")
    call(["rm", "-rf"] + junk)

    # Make sure we use the correct version of Tk:
    tcl_current = framework_dir + "Tcl.framework/Versions/Current"
    tk_current = framework_dir + "Tk.framework/Versions/Current"
    if os.path.exists(tcl_current):
        os.remove(tcl_current)
    os.symlink(tk_ver, tcl_current)
    if os.path.exists(tk_current):
        os.remove(tk_current)
    os.symlink(tk_ver, tk_current)

    # Add a symlink so that py2app 0.13 can find the Tcl Scripts directory
    libdir = "dist/SnapPy.app/Contents/lib/"
    os.mkdir(libdir)
    os.symlink("../Frameworks/Tk.Framework/Versions/Current/Resources/Scripts",
               libdir + "tk%s"%tk_ver)

    # Add a symlink so that Tcl will be able to find its "init.tcl"
    os.symlink('Versions/Current/Resources', framework_dir + 'Tcl.framework/Resources')

def package_app(dmg_name):
    """
    Create a disk image containing the app, with a nice background and
    a symlink to the Applications folder.
    """
    image_dir = "disk_images"
    if not os.path.exists(image_dir):
        os.mkdir(image_dir)
    mount_name = "/Volumes/SnapPy"
    dmg_path = os.path.join(image_dir, dmg_name + ".dmg")
    temp_path = os.path.join(image_dir, dmg_name + "-tmp.dmg")
    # Make sure the dmg isn't currently mounted, or this won't work.  
    while os.path.exists(mount_name):
        print("Trying to eject " + mount_name)
        os.system("hdiutil detach " + mount_name)
    # Remove old dmgs if they exist.
    if os.path.exists(dmg_path):
        os.remove(dmg_path)
    if os.path.exists(temp_path):
        os.remove(temp_path)
    # Add symlink to /Applications if not there.
    if not os.path.exists("dist/Applications"):
        os.symlink("/Applications/", "dist/Applications")

    # Copy over the background and .DS_Store file.
    call(["rm", "-rf", "dist/.background"])
    os.mkdir("dist/.background")
    call(["cp", "dmg-maker/background.png", "dist/.background"])
    call(["cp", "dmg-maker/dotDS_Store", "dist/.DS_Store"])
        
    # Figure out the needed size.
    raw_size, errors = Popen(["du", "-sh", "dist"], stdout=PIPE).communicate()
    size, units = re.search("([0-9.]+)([KMG])", str(raw_size)).groups()
    new_size = "%d" % ceil(1.2 * float(size)) + units
    # Run hdiutil to build the dmg file.:
    call(["hdiutil", "makehybrid", "-hfs", "-hfs-volume-name", "SnapPy",
        "-hfs-openfolder", "dist", "dist", "-o", temp_path])
    call(["hdiutil", "convert", "-format", "UDZO", temp_path, "-o", dmg_path])
    os.remove(temp_path)
    # Delete the symlink to /Applications or egg_info will be glacial on newer setuptools.
    os.remove("dist/Applications")

def do_release(python, dmg_name):
    freshen_SnapPy(python)
    build_app(python)
    cleanup_app(python)
    package_app(dmg_name)

if os.path.exists('/Users/dunfield/pythons'):
    print('Using virtualenv Pythons')
    python2 = '/Users/dunfield/pythons/py27/bin/python'
    python3 = '/Users/dunfield/pythons/py36/bin/python'
else:
    framework = '/Library/Frameworks/Python.framework'
    print('Using python from %s'%framework)
    python2 = os.path.join(framework, 'Versions', '2.7', 'bin', 'python')
    python3 = os.path.join(framework, 'Versions', '3.7', 'bin', 'python3')

do_release(python2, "SnapPy-Python2")
do_release(python3, "SnapPy-Python3")
