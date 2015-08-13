from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os, glob


c_code = glob.glob('c_library/*.c')
cython_code = glob.glob('cython_code/*.pyx')
headers = ['c_headers', 'cython_headers']


ext_modules = []
for cython_file in cython_code:
    ext_modules.append(
        Extension(name = os.path.basename(cython_file)[:-4], 
                  sources = [cython_file] + c_code,
                  include_dirs = headers)
    )
                        
setup(
  ext_modules = cythonize(ext_modules, include_path=headers)
)
