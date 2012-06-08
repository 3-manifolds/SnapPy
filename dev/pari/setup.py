from distutils.core import setup, Command
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os

pari_include_dir = ['.', 'pari-2.5.1/include']
pari_gen = Extension('pari.gen',
                     sources = ['pari/gen.pyx'],
                     include_dirs = pari_include_dir + ['pari', '/usr/local/include'],
                     library_dirs = ['pari-2.5.1/lib'],
                     libraries = ['gmp', 'pari'])

class clean(Command):
    user_options = []
    def initialize_options(self):
        pass 
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -rf build dist *.pyc')
        os.system('rm -rf pari*.egg-info')
        os.system('rm -f */*.c')


setup(
  name = 'pari',
  packages = ['pari'],
  cmdclass = {'build_ext': build_ext, 'clean':clean},
  ext_modules = [pari_gen],
)
