"""
Installation script for the snappy module.

Depends heavily on setuptools.
"""
no_setuptools_message = """
You need to have setuptools installed to build the snappy module, e.g. by:

  curl -O http://peak.telecommunity.com/dist/ez_setup.py
  sudo python ez_setup.py

or by installing the python-setuptools package (Debian/Ubuntu) or
python-setuptools-devel package (Fedora).
"""

no_cython_message = """
You need to have Cython (>= 0.11.2) installed to build the snappy
module, e.g.

  sudo python -m easy_install cython

"""

no_sphinx_message = """
You need to have Sphinx (>= 0.6.1) installed to build the snappy
module, e.g.

  sudo python -m easy_install sphinx

"""

try:
    import setuptools
    import pkg_resources
except ImportError:
    raise ImportError(no_setuptools_message)

# Make sure we have Cython and Sphinx installed before proceeding

try:
    pkg_resources.working_set.require('cython>=0.11.2')
except pkg_resources.DistributionNotFound:
    raise ImportError(no_cython_message)

try:
    pkg_resources.working_set.require('sphinx>=0.6.1')
except pkg_resources.DistributionNotFound:
    raise ImportError(no_sphinx_message)

# Remove '.' from the path so that Sphinx doesn't try to load the SnapPy module directly

import sys, os, glob
try:
    sys.path.remove(os.path.realpath(os.curdir))
except:
    pass

# Build the qd library if necessary

if not os.path.exists('qd') and 'clean' not in sys.argv:
    os.system('cd qd_src ; bash build_qd.bash ')

from distutils.extension import Extension
from Cython.Distutils import build_ext
from setuptools import setup, Command
from pkg_resources import load_entry_point

# A real clean

class clean(Command):
    user_options = []
    def initialize_options(self):
        pass 
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -rf build dist *.pyc')
        os.system('rm -rf snappy*.egg-info')
        os.system('rm -rf snappy/doc')
        os.system('rm */Cy*.c */Cy*.h */Cy*.cpp')
        os.system('rm SnapPy.c SnapPy.h SnapPyHP.cpp SnapPyHP.h')

class build_docs(Command):
    user_options = []
    def initialize_options(self):
        pass 
    def finalize_options(self):
        pass
    def run(self):
        sphinx_cmd = load_entry_point('Sphinx>=0.6.1', 'console_scripts', 'sphinx-build')
        sphinx_args = ['sphinx', '-a', '-E', '-d', 'doc-source/_build/doctrees',
                       'doc-source', 'snappy/doc']
        sphinx_cmd(sphinx_args)

# C source files we provide

base_code = glob.glob(os.path.join('kernel_code','*.c'))
unix_code = glob.glob(os.path.join('unix_kit','*.c'))
unix_code.remove(os.path.join('unix_kit','unix_UI.c'))
unix_code.remove(os.path.join('unix_kit','decode_new_DT.c'))
addl_code = glob.glob(os.path.join('addl_code', '*.c')) + glob.glob(os.path.join('addl_code', '*.cc'))
code  =  base_code + unix_code + addl_code

# Symlinks for building the high precision version

def make_symlinks(source_files, target_dir):
    if sys.platform == 'win32':
        import shutil
        symlink = shutil.copyfile
    else:
        symlink = os.symlink
    if not os.path.exists(target_dir):
        os.mkdir(target_dir)
    for path in source_files:
        filename, ext = os.path.splitext(os.path.basename(path))
        new_ext = '.cpp' if ext == '.c' else ext
        link_name = os.path.join(target_dir, filename + new_ext)
        if sys.platform == 'win32':
            source = path
        else:
            source = os.path.join('..', path)
        if not os.path.exists(link_name):
            print 'linking %s -> %s'%(link_name, source) 
            symlink(source, link_name)

def setup_symlinks():
    make_symlinks(base_code, 'hp_kernel_code')
    make_symlinks(unix_code, 'hp_unix_kit')                                 
    unix_headers = glob.glob(os.path.join('unix_kit', '*.h'))
    make_symlinks(unix_headers, 'hp_unix_kit')
    make_symlinks(addl_code, 'hp_addl_code')
    addl_headers = glob.glob(os.path.join('addl_code', '*.h'))
    make_symlinks(addl_headers, 'hp_addl_code')
    headers = glob.glob(os.path.join('headers', '*.h'))
    make_symlinks(headers, 'hp_headers')                                 

class build_symlinks(Command):
    user_options = []
    def initialize_options(self):
        pass 
    def finalize_options(self):
        pass
    def run(self):
        setup_symlinks()

# C++ source files we provide

if not os.path.exists('hp_kernel_code') and 'clean' not in sys.argv:
    setup_symlinks()

hp_base_code = glob.glob(os.path.join('hp_kernel_code','*.cpp'))
hp_unix_code = glob.glob(os.path.join('hp_unix_kit','*.cpp'))
hp_addl_code = glob.glob(os.path.join('hp_addl_code', '*.cpp'))# + glob.glob(os.path.join('addl_code', '*.cc'))
hp_code  =  hp_base_code + hp_unix_code + hp_addl_code

# We replace the SnapPea kernel module Dirichlet_precision.c,
# so let's not link against it.
#code.remove(os.path.join('kernel_code','Dirichlet_precision.c'))

# The SnapPy extension
SnapPyC = Extension(
    name = 'snappy.SnapPy',
    sources = ['SnapPy.pyx', 'SnapPycore.pxi', 'SnapPy.pxi'] + code, 
    include_dirs = ['headers', 'unix_kit', 'addl_code'],
    extra_objects = [])

# The high precision SnapPy extension
SnapPyHP = Extension(
    name = 'snappy.SnapPyHP',
    sources = ['SnapPyHP.pyx', 'SnapPycore.pxi', 'SnapPy.pxi'] + hp_code, 
    include_dirs = ['hp_headers', 'hp_unix_kit', 'hp_addl_code', 'qd/include/'],
    libraries = ['qd'],
    library_dirs = ['qd/lib'],
    language='c++',
    extra_objects = [])

# The CyOpenGL extension
CyOpenGL_includes = ['.']
CyOpenGL_libs = []
CyOpenGL_extras = []
CyOpenGL_extra_link_args = []
if sys.platform == 'darwin':
    CyOpenGL_includes += ['/System/Library/Frameworks/OpenGL.framework/Versions/Current/Headers/']
    CyOpenGL_extra_link_args = ['-framework', 'OpenGL']
elif sys.platform == 'linux2':
    CyOpenGL_includes += ['/usr/include/GL']
    CyOpenGL_libs += ['GL', 'GLU']
elif sys.platform == 'win32':
    CyOpenGL_includes += ['/mingw/include/GL']
    CyOpenGL_extras += ['/mingw/lib/libopengl32.a',
                        '/mingw/lib/libglu32.a']
    

CyOpenGL = Extension(
    name = 'snappy.CyOpenGL',
    sources = ['opengl/CyOpenGL.pyx'], 
    include_dirs = CyOpenGL_includes,
    libraries = CyOpenGL_libs,
    extra_objects = CyOpenGL_extras,
    extra_link_args = CyOpenGL_extra_link_args)

# Twister

twister_main_path = 'Twister/lib/'
twister_main_src = [twister_main_path + 'py_wrapper.cpp']
twister_kernel_path = twister_main_path + 'kernel/'
twister_kernel_src = [twister_kernel_path + file for file in
                      ['twister.cpp', 'manifold.cpp', 'parsing.cpp', 'global.cpp']]

TwisterCore = Extension(
	name = 'snappy.twister.twister_core',
	sources = twister_main_src + twister_kernel_src,
	include_dirs=[twister_kernel_path],
	language='c++' )

#ext_modules = [SnapPyC, CyOpenGL, TwisterCore]
#ext_modules = [SnapPyHP, CyOpenGL, TwisterCore]
ext_modules = [SnapPyC, SnapPyHP, CyOpenGL, TwisterCore]

try:
    import sage
    install_requires = ['plink>=1.6', 'ipython', 'pypng', 'spherogram>=1.3']
except ImportError:
    install_requires = ['plink>=1.6', 'ipython>=0.13', 'pypng', 'spherogram>=1.3', 'pyttk', 'cypari>=1.0']
    if sys.platform == 'win32':
        install_requires.append('pyreadline>=2.0')
    
    
# Get version number:

exec(open('snappy/version.py').read())


# Off we go ...
setup( name = 'snappy',
       version = version,
       zip_safe = False,
       install_requires = install_requires,
       dependency_links = ['http://www.math.uic.edu/t3m/plink/',
                           'http://www.math.uic.edu/t3m/SnapPy-nest'],
       packages = ['snappy', 'snappy/manifolds', 'snappy/twister',
                   'snappy/snap', 'snappy/snap/t3mlite', 'snappy/ptolemy'],
       package_data = {
        'snappy' : ['togl/*-tk*/Togl2.0/*',
                    'togl/*-tk*/Togl2.1/*',
                    'togl/*-tk*/mactoolbar*/*',
                    'info_icon.gif', 'SnapPy.ico',
                    'doc/*.*',
                    'doc/_images/*',
                    'doc/_sources/*',
                    'doc/_static/*'],
        'snappy/manifolds' : ['manifolds.sqlite',
                              'more_manifolds.sqlite',
                              'HTWKnots/*.gz'],
        'snappy/twister' : ['surfaces/*'],
        'snappy/ptolemy':['testing_files/*magma_out.bz2', 'testing_files_generalized/*magma_out.bz2'],
        },
       package_dir = {'snappy/twister':'Twister/lib'},
       ext_modules = ext_modules,
       cmdclass =  {'build_ext': build_ext,
                    'clean' : clean,
                    'build_symlinks': build_symlinks,
                    'build_docs': build_docs},
       entry_points = {'console_scripts': ['SnapPy = snappy.app:main']},
       author = 'Marc Culler and Nathan Dunfield',
       author_email = 'culler@math.uic.edu, nmd@illinois.edu',
       description = "Python application based on Jeff Weeks' SnapPea",
       license = 'GPL v2 or later',
       keywords = 'hyperbolic 3-manifolds',
       url = 'http://www.math.uic.edu/t3m',
       )
