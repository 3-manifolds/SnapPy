from distutils.core import setup
import py2exe
import glob
import IPython
import os

APP = [{
        'script' : 'SnapPy.py',
        'icon_resources' : [(1, 'SnapPy.ico')]
       }]

DATA_FILES = [
('snappy', ['../snappy/info_icon.gif', '../snappy/SnapPy.ico']), 
('snappy/manifolds', ['../snappy/manifolds/manifolds.sqlite',
                      '../snappy/manifolds/more_manifolds.sqlite']),
('snappy/manifolds/HTWKnots',
 glob.glob('../snappy/manifolds/HTWKnots/*.*')),
('snappy/togl/win32-tk8.5', ['../snappy/togl/win32-tk8.5/Toglstub20.lib']),
('snappy/togl/win32-tk8.5/Togl2.0',
glob.glob('../snappy/togl/win32-tk8.5/Togl2.0/*')),
('snappy/twister/surfaces',
 glob.glob('../Twister/lib/surfaces/*')),
('snappy/doc',
glob.glob('../snappy/doc/*.*')),
('snappy/doc/_static',
glob.glob('../snappy/doc/_static/*.*')),
('snappy/doc/_images',
glob.glob('../snappy/doc/_images/*.*')),
('snappy/doc/_sources',
glob.glob('../snappy/doc/_sources/*.*')),
('IPython/config/profile',
[os.path.join(IPython.__path__[0], 'config', 'profile', 'README_STARTUP')] )
]

OPTIONS = {
'excludes':
  'scipy,numpy',
'packages': 
  'snappy,snappy.manifolds,snappy.twister,IPython,plink,readline,pyreadline',
'includes':
  'gzip,tarfile,pydoc,pyx',
'skip_archive':
  1,
'dist_dir':
  'SnapPy',
}

setup(
    windows=APP,
    data_files=DATA_FILES,
    options={'py2exe': OPTIONS},
)
