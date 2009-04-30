import SnapPy
#import os, sys, operator, types, re, gzip, struct, tempfile, tarfile
#import signal
#import plink
#import SnapPy.polyviewer
#import SnapPy.horoviewer
#import SnapPy.manifolds
#import ctypes
#import ctypes.util

m = SnapPy.Manifold("m004")
d = SnapPy.DirichletDomain(m)
d.view()
d.viewer.window.mainloop()

#__requires__ = 'ipython==0.9.1'
#import pkg_resources
#from pkg_resources import load_entry_point
#
#sys.exit(
#   load_entry_point('ipython==0.9.1', 'console_scripts', 'ipython')()
#)
