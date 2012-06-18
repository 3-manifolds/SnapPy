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

# Make sure we have Cython installed before proceeding

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
sys.path.remove(os.path.realpath(os.curdir))

# Hack to patch setuptools so that it treats Cython
# as a replacement for pyrex.

from distutils.core import Extension as _Extension
from setuptools.dist import _get_unpatched
_Extension = _get_unpatched(_Extension)

try:
    from Cython.Distutils import build_ext
except ImportError:
    have_cython = False
else:
    have_cython = True


class Extension(_Extension):
    """
    This modified version of setuptools Extension allows us
    to use Cython instead of pyrex.  If Cython is not installed
    on this system, it will assume that a Cython-generated .c
    file is present in the distribution.
    """
    if not have_cython:
        # convert .pyx extensions to .c
        def __init__(self,*args,**kw):
            _Extension.__init__(self,*args,**kw)
            sources = []
            for s in self.sources:
                if s.endswith('.pyx'):
                    sources.append(s[:-3]+'c')
                else:
                    sources.append(s)
            self.sources = sources

class Library(Extension):
    """Just like a regular Extension, but built as a library instead"""

import sys, distutils.core, distutils.extension
distutils.core.Extension = Extension
distutils.extension.Extension = Extension
if 'distutils.command.build_ext' in sys.modules:
    sys.modules['distutils.command.build_ext'].Extension = Extension
# End of hack

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
        os.system('rm */Cy*.c */Cy*.h */Cy*.cpp SnapPy.c SnapPy.h')

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
addl_code = glob.glob(os.path.join('addl_code', '*.c')) + glob.glob(os.path.join('addl_code', '*.cc'))
code  =  base_code + unix_code + addl_code

# We replace the SnapPea kernel module Dirichlet_precision.c,
# so let's not link against it.
code.remove(os.path.join('kernel_code','Dirichlet_precision.c'))

# The SnapPy extension
SnapPyC = Extension(
    name = 'snappy.SnapPy',
    sources = ['SnapPy.pxi','SnapPy.pyx'] + code, 
    include_dirs = ['headers', 'unix_kit', 'addl_code'],
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

ext_modules = [SnapPyC, CyOpenGL, TwisterCore]

try:
    import sage
    install_requires = ['plink>=1.2', 'ipython', 'pypng']
except ImportError:
    install_requires = ['plink>=1.2', 'ipython>=0.12', 'pypng', 'pyttk']
    
# Get version number:

exec(open('snappy/version.py').read())


# Off we go ...
setup( name = 'snappy',
       version = version,
       zip_safe = False,
       install_requires = install_requires,
       dependency_links = ['http://math.uic.edu/t3m/plink/', 'http://math.uic.edu/t3m/SnapPy/'],
       packages = ['snappy', 'snappy/manifolds', 'snappy/twister'],
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
                              'HTWKnots/*.gz',
                              'MTLinks/*.gz'],
        'snappy/twister' : ['surfaces/*'], 
        },
       package_dir = {'snappy/twister':'Twister/lib'},
       ext_modules = ext_modules,
       cmdclass =  {'build_ext': build_ext, 'clean' : clean, 'build_docs': build_docs},
       entry_points = {'console_scripts': ['SnapPy = snappy.app:main']},
       author = 'Marc Culler and Nathan Dunfield',
       author_email = 'culler@math.uic.edu, nmd@illinois.edu',
       description = "Python application based on Jeff Weeks' SnapPea",
       license = 'GPL',
       keywords = 'hyperbolic 3-manifolds',
       url = 'http://www.math.uic.edu/t3m',
       )
