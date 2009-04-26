from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from distutils.command.install_data import install_data
import os, glob

# The default is to build pari inside this directory,
# but you can modify this either here or by creating
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


class SnapPy_install_data(install_data):
    def finalize_options (self):
        self.set_undefined_options('install',
                                   ('install_lib', 'install_dir'),
                                   ('root', 'root'),
                                   ('force', 'force'),
                                   )
    
base_code = glob.glob(os.path.join("kernel_code","*.c"))
unix_code = glob.glob(os.path.join("unix_kit","*.c"))
unix_code.remove(os.path.join("unix_kit","unix_UI.c"))
addl_code = glob.glob(os.path.join("addl_code", "*.c")) + glob.glob(os.path.join("addl_code", "*.cc"))
code  =  base_code + unix_code + addl_code

data_dir = "SnapPy/manifolds"
closed = glob.glob(os.path.join(data_dir,"ClosedCensusData","Cl*"))
cusped = glob.glob(os.path.join(data_dir,"CuspedCensusData","t*"))
knots  = glob.glob(os.path.join(data_dir,"HTWKnots","*.gz"))
links  = [os.path.join(data_dir,"ChristyLinks.tgz")]

togl_dirs = glob.glob(os.path.join("SnapPy","*-tk*","Togl2.0"))
togl_files = []
for dir in togl_dirs:
    filelist = glob.glob(os.path.join(dir,"*"))
    togl_files.append((dir, filelist))

SnapPyC = Extension(name = "SnapPy.SnapPy",
                   sources = ["SnapPy.pxi","SnapPy.pyx"] + code, 
                   include_dirs = ["headers", "unix_kit", "addl_code"] + pari_include_dir,
                   extra_objects = [] + pari_extra_objects)

setup( name = "SnapPy",
       ext_modules = [SnapPyC],
       packages = ["SnapPy", "SnapPy/manifolds"],
       cmdclass = {'build_ext': build_ext, "install_data" : SnapPy_install_data},
       data_files = [(data_dir+"/ClosedCensusData", closed),
                     (data_dir+"/CuspedCensusData", cusped),
                     (data_dir+"/HTWKnots", knots),
                     (data_dir, links)] + togl_files
       )



       


