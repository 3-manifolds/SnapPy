from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

pari_gen = Extension('pari.gen',
              sources = ['src/gen.pyx'],
              libraries = ['pari', 'gmp'])
setup(
  name = 'pari',
  packages = ['pari'],
  package_dir = {'pari':'src'},
  cmdclass = {'build_ext': build_ext},
  ext_modules = [pari_gen],
)
