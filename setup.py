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
    raise ImportError, no_setuptools_message

# Make sure we have Cython installed before proceeding

try:
    pkg_resources.working_set.require("cython>=0.11.2")
except pkg_resources.DistributionNotFound:
    raise ImportError, no_cython_message

try:
    pkg_resources.working_set.require("sphinx>=0.6.1")
except pkg_resources.DistributionNotFound:
    raise ImportError, no_sphinx_message

# Remove "." from the path so that Sphinx doesn't try to load the SnapPy module directly

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
        os.system("rm -rf build dist *.pyc")
        os.system("rm -rf snappy*.egg-info")
        os.system("rm -rf snappy/doc")
        for file in glob.glob("*.pyx"):
            os.system("rm -rf " + file[:-3] + "c")

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

# The default is to build pari inside this directory,
# but you can modify this either here, or by creating
# a file pari_path which overides them.  

pari_include_dir = ["pari-2.3.4/include/", "pari-2.3.4/include/pari"]
pari_extra_objects = ["pari-2.3.4/lib/libpari.a"]

try:
    from pari_paths import *
except:
    pass

# C source files we provide
base_code = glob.glob(os.path.join("kernel_code","*.c"))
unix_code = glob.glob(os.path.join("unix_kit","*.c"))
unix_code.remove(os.path.join("unix_kit","unix_UI.c"))
addl_code = glob.glob(os.path.join("addl_code", "*.c")) + glob.glob(os.path.join("addl_code", "*.cc"))
code  =  base_code + unix_code + addl_code

try:
    import sage
    snappy_extra_objects = []
except ImportError:
    snappy_extra_objects = pari_extra_objects

# The SnapPy extension
SnapPyC = Extension(
    name = "snappy.SnapPy",
    sources = ["SnapPy.pxi","SnapPy.pyx"] + code, 
    include_dirs = ["headers", "unix_kit", "addl_code"] + pari_include_dir,
    extra_objects = snappy_extra_objects)

# The CyOpenGL extension
CyOpenGL_includes = []
CyOpenGL_libs = []
CyOpenGL_extras = []
CyOpenGL_extra_link_args = []
if sys.platform == 'darwin':
    CyOpenGL_includes += ['/System/Library/Frameworks/OpenGL.framework/Versions/Current/Headers/']
    CyOpenGL_extra_link_args = ["-framework", "OpenGL"]
elif sys.platform == 'linux2':
    CyOpenGL_includes += ['/usr/include/GL']
    CyOpenGL_libs += ['GL', 'GLU']
elif sys.platform == 'win32':
    CyOpenGL_includes += ['/mingw/include/GL']
    CyOpenGL_extras += ['/mingw/lib/libopengl32.a',
                        '/mingw/lib/libglu32.a']
    

CyOpenGL = Extension(
    name = "snappy.CyOpenGL",
    sources = ["CyOpenGL.pyx"], 
    include_dirs = CyOpenGL_includes,
    libraries = CyOpenGL_libs,
    extra_objects = CyOpenGL_extras,
    extra_link_args = CyOpenGL_extra_link_args)

# We use PARI to augment the SnapPea kernels ability to compute
# homology.  If we're within Sage, we use that instead.

CyPari = Extension(
    name = "snappy.CyPari",
    sources = ["CyPari.pyx"], 
    include_dirs = pari_include_dir, 
    extra_objects = pari_extra_objects,
)

# Twister

twister_dir = "Twister"
CyTwister = Extension(
    name = "snappy.CyTwister",
    language="c++",
    sources = ["CyTwister.pyx"] + glob.glob(twister_dir + "/*.cpp"),
    include_dirs = [twister_dir]
)

try:
    import sage
    ext_modules = [SnapPyC, CyOpenGL]
except ImportError:
    ext_modules = [SnapPyC, CyOpenGL, CyPari, CyTwister]

# Get version number:

execfile('snappy/version.py')

# Off we go ...
setup( name = "snappy",
       version = version,
       zip_safe = False,
       install_requires = ['plink>=1.1', 'ipython>=0.9'],
       dependency_links = ['http://math.uic.edu/t3m/plink/', 'http://math.uic.edu/t3m/SnapPy/'],
       packages = ["snappy", "snappy/manifolds"],
       package_data = {
        'snappy' : ['*-tk*/Togl2.0/*',
                    'doc/*.*',
                    'doc/_images/*',
                    'doc/_sources/*',
                    'doc/_static/*'],
        'snappy/manifolds' : ['CensusKnots.tgz',
                              'ChristyLinks.tgz',
                              'morwen8.tgz',
                              'ClosedCensusData/*.txt',
                              'CuspedCensusData/*.bin',
                              'HTWKnots/*.gz',
                              'MTLinks/*.gz']
        },
       ext_modules = ext_modules,
       cmdclass =  {'build_ext': build_ext, 'clean' : clean, 'build_docs': build_docs},
       entry_points = {'console_scripts': ['SnapPy = snappy.app:main']},
       author = "Marc Culler and Nathan Dunfield",
       author_email = "culler@math.uic.edu, nmd@illinois.edu",
       description = "Python application based on Jeff Weeks' SnapPea",
       license = "GPL",
       keywords = "hyperbolic 3-manifolds",
       url = "http://www.math.uic.edu/t3m",
       )

