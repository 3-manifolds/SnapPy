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

from setuptools import setup
import os, glob

# The default is to build pari inside this directory,
# but you can modify this either here, or by creating
# a file pari_path which overides them.  

pari_include_dir = ["pari-2.3.4/include/", "pari-2.3.4/include/pari"]
pari_extra_objects = ["pari-2.3.4/lib/libpari.a"]

# If we're being called from SAGE, we just want to use it's copy of PARI
try:
    import sage
    sage_root = os.environ["SAGE_ROOT"]
    pari_include_dir = [sage_root + "/local/include/pari"]
    pari_extra_objects = [sage_root + "/local/lib/libpari.a",] + glob.glob(sage_root +  "/pkgs/sage/local/lib/libgmp.*")
except:
    pass

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

# The SnapPy extension
SnapPyC = Extension(
    name = "snappy.SnapPy",
    sources = ["SnapPy.pxi","SnapPy.pyx"] + code, 
    include_dirs = ["headers", "unix_kit", "addl_code"] + pari_include_dir,
    extra_objects = [] + pari_extra_objects)

# Off we go ...
setup( name = "snappy",
       version = "1.0a",
       zip_safe = False,
       install_requires = ['ipython>=0.9', 'PyOpenGL>2.9'],
       packages = ["snappy", "snappy/manifolds"],
       package_data = {
        'snappy' : ['*-tk*/Togl2.0/*',
                    'doc/*'],
        'snappy/manifolds' : ['ChristyLinks.tgz',
                              'ClosedCensusData/*.txt',
                              'CuspedCensusData/*.bin',
                              'HTWKnots/*.gz']
        },
       ext_modules = [SnapPyC],
       cmdclass =  {'build_ext': build_ext},
       author = "Marc Culler and Nathan Dunfield",
       author_email = "culler@math.uic.edu, nmd@illinois.edu",
       description = "Python application based on Jeff Weeks' SnapPea",
       license = "GPL",
       keywords = "hyperbolic 3-manifolds",
       url = "http://www.math.uic.edu/~t3m",
       download_url = ""
       )

