from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os, glob

code_dirs = [ ('.c', 'c_library'), ('.pyx', 'cython_code'), ('.pxd','cython_code')]

sources = []
for  code_type, directory in code_dirs:
    sources += glob.glob(os.path.join(directory, '*'+code_type))

print sources




ext_modules = [Extension(name = "cytest",
                         sources = sources,
                         include_dirs = ['headers']
                         )]

setup(
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
