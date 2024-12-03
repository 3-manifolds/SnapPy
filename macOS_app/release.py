#! /usr/bin/env python
import os
import sys
import re
import shutil
from glob import glob
import subprocess
import configparser
from subprocess import check_call, call, Popen, PIPE
from math import ceil
from check_target import TkChecker
from notabot import Notarizer

PYTHON_VERSION = '3.13'
PYTHON_VERSION_SHORT = PYTHON_VERSION.replace('.', '')
APP_PYTHON = 'python' + PYTHON_VERSION
PYTHON_ZIP = 'python' + PYTHON_VERSION_SHORT + '.zip'
FRAMEWORKS_TARBALL = 'Frameworks-' + PYTHON_VERSION + '.tgz'

# Make sure that we have our frameworks.
if not os.path.exists(FRAMEWORKS_TARBALL):
    print("Please build the frameworks for SnapPy.app first.")
    sys.exit(1)

os.environ['_PYTHON_HOST_PLATFORM'] = 'macosx-10.13-universal2'
os.environ['ARCHFLAGS'] = '-arch arm64 -arch x86_64'

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
    os.chdir("macOS_app")
    
def build_app(python):
    """
    Build the standalone app bundle.
    """
    check_call([python, "setup.py", "py2app"])
    # Replace the frameworks that py2app installs with our own signed frameworks.
    shutil.rmtree(os.path.join('dist', 'SnapPy.app', 'Contents', 'Frameworks'))
    check_call(['tar', 'xfz', FRAMEWORKS_TARBALL])
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
    contents = os.path.join('dist', 'SnapPy.app', 'Contents')
    frameworks = os.path.join(contents, 'Frameworks')
    resources = os.path.join(contents, 'Resources')
    python_lib_dir = os.path.join(frameworks, 'Python.framework', 'Versions',
        'Current', 'lib', python)
    # Remove the dev directory
    dev_directory = os.path.join(resources, 'lib', python, 'snappy', 'dev')
    shutil.rmtree(dev_directory, ignore_errors=True)
    # Remove Python.org junk that may or may not have been added by py2app.
    shutil.rmtree(os.path.join(resources, 'lib', 'tcl8.6'), ignore_errors=True)
    shutil.rmtree(os.path.join(resources, 'lib', 'tcl8'), ignore_errors=True)
    shutil.rmtree(os.path.join(resources, 'lib', 'tk8.6'), ignore_errors=True)
    # Remove everything in the python lib except for the lib-dynloads directory,
    # since all of that stuff is in the zip filesystem (as bytecode).
    tmp_dir = os.path.join(frameworks, 'Python.framework', 'Versions', 'Current',
        'lib', 'tmp')
    os.mkdir(tmp_dir)
    shutil.move(os.path.join(python_lib_dir, 'lib-dynload'), tmp_dir)
    shutil.rmtree(python_lib_dir, ignore_errors= True)
    os.rename(tmp_dir, python_lib_dir)

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

def do_release(python, dmg_name, freshen=True):
    if freshen:
        freshen_SnapPy(python)
    build_app(python)
    cleanup_app(python)

class SnapPyNotarizer(Notarizer):

    def build_dmg(self):
        app_name = self.config['app']['app_name']
        dmg_file = self.config['app']['dmg_path']
        print('building dmg %s for %s'%(dmg_file, app_name))
        package_app(app_name)

if __name__ == '__main__':
    freshen = '--no-freshen' not in sys.argv
    do_release(APP_PYTHON, "SnapPy", freshen)
    if '--notarize' in sys.argv:
        notarizer = SnapPyNotarizer('notabot.cfg')
        notarizer.run()
        
