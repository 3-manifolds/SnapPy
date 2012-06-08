from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

pari_include_dir = ['.', 'pari-2.5.1/include']
pari_extra_objects = ['pari-2.5.1/lib/libpari.a']

pari_gen = Extension('pari.gen',
              sources = ['src/gen.pyx'],
              include_dirs = pari_include_dir, 
              extra_objects = pari_extra_objects,
              libraries = ['gmp'])
setup(
  name = 'pari',
  packages = ['pari'],
  package_dir = {'pari':'src'},
  cmdclass = {'build_ext': build_ext},
  ext_modules = [pari_gen],
)
