from distutils.core import setup
import py2exe
import glob

APP = [{
        'script' : 'SnapPy.py',
        'icon_resources' : [(1, 'SnapPy.ico')]
       }]

DATA_FILES = [
('snappy/manifolds',
 glob.glob('../snappy/manifolds/*.tgz')),
('snappy/manifolds/HTWKnots',
 glob.glob('../snappy/manifolds/HTWKnots/*.*')),
('snappy/manifolds/ClosedCensusData',
 glob.glob('../snappy/manifolds/ClosedCensusData/*.*')),
('snappy/manifolds/CuspedCensusData',
 glob.glob('../snappy/manifolds/CuspedCensusData/*.*')),
('snappy/win32-tk8.5/Togl2.0',
glob.glob('../snappy/win32-tk8.5/Togl2.0/*')),
('OpenGL/DLLS',
['glut32.dll']),
('snappy/doc',
glob.glob('../snappy/doc/*.*')),
('snappy/doc/_static',
glob.glob('../snappy/doc/_static/*.*')),
('snappy/doc/_images',
glob.glob('../snappy/doc/_images/*.*')),
('snappy/doc/_sources',
glob.glob('../snappy/doc/_sources/*.*')),
]

OPTIONS = {
'excludes':
  'scipy,numpy',
'packages': 
  'snappy,snappy.manifolds,IPython,plink',
'includes':
  'gzip,tarfile,readline,pydoc',
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
