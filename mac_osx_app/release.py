#! /usr/bin/env python
import os
import sys
import re
import shutil
from glob import glob
from subprocess import check_call, Popen, PIPE
from math import ceil
from check_target import TkChecker

try:
    import pyx
except ImportError:
    print("ERROR: Need to install PyX!")
    sys.exit(1)

def get_tk_ver(python):
    """
    Figure out which version of Tk is used by this python.
    """
    out, errors = Popen([python, "-c", "import _tkinter; print(_tkinter.TK_VERSION)"],
                        stdout=PIPE, text=True).communicate()
    return out.strip()

def freshen_SnapPy(python):
    """
    Build SnapPy and install it twice to make sure the documentation is
    up to date.
    """
    os.chdir("../")
    check_call(["git", "pull"])
    check_call([python, "setup.py", "pip_install"])
    os.chdir("mac_osx_app")
    
def build_app(python):
    """
    Build the standalone app bundle.
    """
    try:
        checker = TkChecker()
        if checker.Tk_target != '10.9' or checker.Tcl_target != '10.9':
            print(checker)
            raise RuntimeError('Tk was not built for macOSX 10.9!')
        else:
            print('Tk looks fine on this build system.')
    except:
        raise RuntimeError('Tk was not built for macOSX 10.9!')
    check_call([python, "setup.py", "py2app"])
    # Replace the frameworks that py2app installs with our own signed frameworks.
    shutil.rmtree(os.path.join('dist', 'SnapPy.app', 'Contents', 'Frameworks'))
    check_call(['tar', 'xfz', 'Frameworks.tgz'])
    contents = os.path.join('dist', 'SnapPy.app', 'Contents')
    resources = os.path.join(contents, 'Resources')
    frameworks = os.path.join(contents, 'Frameworks')
    os.rename('Frameworks', frameworks)
    shutil.copy('Info.plist', contents)
    shutil.copy('__boot__.py', resources)

def cleanup_app(python):
    """
    Tidy things up.
    """
    extra_dynload = glob('dist/SnapPy.app/Contents/Resources/lib/python*/lib-dynload')[0]
    shutil.rmtree(extra_dynload)

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
    check_call(["rm", "-rf", "dist/.background"])
    os.mkdir("dist/.background")
    check_call(["cp", "dmg-maker/background.png", "dist/.background"])
    check_call(["cp", "dmg-maker/dotDS_Store", "dist/.DS_Store"])
        
    # Figure out the needed size.
    raw_size, errors = Popen(["du", "-sh", "dist"], stdout=PIPE).communicate()
    size, units = re.search("([0-9.]+)([KMG])", str(raw_size)).groups()
    new_size = "%d" % ceil(1.2 * float(size)) + units
    # Run hdiutil to build the dmg file.:
    check_call(["hdiutil", "makehybrid", "-hfs", "-hfs-volume-name", "SnapPy",
                "-hfs-openfolder", "dist", "dist", "-o", temp_path])
    check_call(["hdiutil", "convert", "-format", "UDZO", temp_path, "-o", dmg_path])
    os.remove(temp_path)
    # Delete the symlink to /Applications or egg_info will be glacial on newer setuptools.
    os.remove("dist/Applications")

def do_release(python, dmg_name):
    freshen_SnapPy(python)
    build_app(python)
    cleanup_app(python)
#    package_app(dmg_name)

if '-m' in sys.argv or '--manual' in sys.argv:
    do_release(sys.executable, "SnapPy-Python" + repr(sys.version_info[0]))
else:
    nmd_python_dir = os.environ['HOME'] + '/pkgs/pythons'
    if os.path.exists(nmd_python_dir):
        print('Using virtualenv Pythons')
        python3 = nmd_python_dir + '/py39/bin/python'
    else:
        python3 = 'python3'
    do_release(python3, "SnapPy")

