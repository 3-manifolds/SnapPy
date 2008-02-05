from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os, glob
    
base_code = glob.glob(os.path.join("kernel_code","*.c"))
unix_code = glob.glob(os.path.join("unix_kit","*.c"))
addl_code = glob.glob(os.path.join("addl_code", "*.c"))
code  =  base_code + unix_code + addl_code
#data_dir ="../SnapPea/manifolds"
#links  = glob.glob(os.path.join("../SnapPea","manifolds","ChristyLinks","L*"))
#closed = glob.glob(os.path.join("../SnapPea","manifolds","ClosedCensusData","Cl*"))
#cusped = glob.glob(os.path.join("../SnapPea","manifolds","CuspedCensusData","t*"))
#knots  = glob.glob(os.path.join("../SnapPea","manifolds","HTWKnots","*.gz"))

# for debugging add 'extra_compile_args = ["-g"]' to the below 

SnapPy = Extension(name = "SnapPy",
                   sources = ["SnapPy.pyx"] + code, 
                   include_dirs = ["headers", "unix_kit", "pari/include/pari"],
                   extra_objects = ['pari/lib/libpari.a'])

setup( name = "SnapPy",
       ext_modules = [SnapPy],
       cmdclass = {'build_ext': build_ext}
       )



       


