from setuptools import setup, Command
import os, sys, glob, plistlib
from py2app import __version__ as py2app_version
import py2app.util

def codesign_adhoc(bundle):
  print('Monkey patching py2app to prevent codesigning, since we do this ourselves at the end')
py2app.util.codesign_adhoc = codesign_adhoc
  
from snappy.version import version as SnapPy_version
try:
  import pyx
except ImportError:
  raise ValueError('Please install PyX.')

class clean(Command):
    user_options = []
    def initialize_options(self):
        pass 
    def finalize_options(self):
        pass
    def run(self):
        os.system("rm -rf build dist *.pyc")
        
if sys.version_info >= (3,4):
    with open('Info.plist', 'rb') as info_plist:
        plist_dict = plistlib.load(info_plist)
else:
    plist_dict = plistlib.readPlist('Info.plist')

runtime_path = os.path.join('@executable_path', os.path.pardir, 'Frameworks', 'Python.framework',
                            'Versions', '%s.%s'%sys.version_info[:2], 'Python')
plist_dict['PyRuntimeLocations'] = [runtime_path]
plist_dict['CFBundleShortVersionString'] = SnapPy_version
plist_dict['CFBundleVersion'] = SnapPy_version
plist_dict['PythonInfoDict']['PythonExecutable'] = sys.executable
plist_dict['PythonInfoDict']['PythonLongVersion'] = sys.version
plist_dict['PythonInfoDict']['py2app']['version'] = py2app_version
packages = 'snappy,spherogram,snappy_manifolds,snappy_15_knots,plink,cypari,FXrays,'
packages += 'IPython,pygments,plink,pyx,knot_floer_homology,low_index,tkinter_gl'

try:
  import jedi
  packages += ',jedi'
except ImportError:
  pass

try:
  import parso
  packages += ',parso'
except ImportError:
  pass

APP = ['SnapPyApp.py']
DATA_FILES = ['SnapPy.sdef']
OPTIONS = {'argv_emulation': False,
           'excludes': 'scipy,numpy,wx,wxversion,wxPython,matplotlib,sphinx,'
                       'idlelib,docutils,curses,cython,Cython,pandas,OpenGL,'
                       'setuptools,test',
           'packages': packages,
           'includes': 'tkinter,gzip,tarfile,readline,pydoc,fractions,pickleshare',
           'iconfile': 'icons/SnapPy.icns',
           'plist'   : plist_dict,
           'arch'    : 'universal2'
}

# Setting custom commands now causes setup.py py2app to crash.
if 'clean' in sys.argv:
    setup(
        app=APP,
        data_files=DATA_FILES,
        cmdclass={'clean': clean}
        )
else:
    setup(
        app=APP,
        data_files=DATA_FILES,
        options={'py2app': OPTIONS},
        setup_requires=['py2app']
        )
