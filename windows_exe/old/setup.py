from distutils.core import setup
import py2exe
import glob
import IPython
import os
import pyx

APP = [{
        'script' : 'SnapPy.py',
        'icon_resources' : [(1, 'SnapPy.ico')]
       }]

DATA_FILES = [
('snappy', ['../python/info_icon.gif', '../python/SnapPy.ico']), 
('snappy/manifolds', ['../python/manifolds/manifolds.sqlite',
                      '../python/manifolds/more_manifolds.sqlite']),
('snappy/manifolds/HTWKnots',
 glob.glob('../python/manifolds/HTWKnots/*.*')),
('snappy/togl/win32VC-tk8.5', ['../python/togl/win32VC-tk8.5/Toglstub21.lib']),
('snappy/togl/win32VC-tk8.5/Togl2.1',
glob.glob('../python/togl/win32VC-tk8.5/Togl2.1/*')),
('snappy/twister/surfaces',
 glob.glob('../Twister/lib/surfaces/*')),
('snappy/doc',
glob.glob('../python/doc/*.*')),
('snappy/doc/_static',
glob.glob('../python/doc/_static/*.*')),
('snappy/doc/_images',
glob.glob('../python/doc/_images/*.*')),
('snappy/doc/_sources',
glob.glob('../python/doc/_sources/*.*')),
('IPython/config/profile',
[os.path.join(IPython.__path__[0], 'config', 'profile', 'README_STARTUP')],
 ),
('IPython/html/static/custom',
glob.glob(os.path.join(IPython.__path__[0], 'html', 'static', 'custom', '*.*'))),
('pyx/data', [os.path.dirname(pyx.__file__) + '/data/pyxrc']), 
('pyx/data/afm', glob.glob(os.path.dirname(pyx.__file__) + '/data/afm/*')), 
('pyx/data/def', glob.glob(os.path.dirname(pyx.__file__) + '/data/def/*')), 
('pyx/data/lfs', glob.glob(os.path.dirname(pyx.__file__) + '/data/lfs/*')), 
]

OPTIONS = {
    'excludes': 'scipy,numpy,iPython.config',
    'packages': 'snappy,snappy.manifolds,snappy.twister,IPython,IPython.html,plink,readline,pyreadline,pyx,lib2to3,cypari.version,spherogram.version',
    'includes': 'gzip,tarfile,pydoc',
    'skip_archive': 1,
    'dist_dir': 'SnapPy',
}

setup(
    windows=APP,
    data_files=DATA_FILES,
    options={'py2exe': OPTIONS},
)
