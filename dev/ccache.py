"""
Speeding up compilation of SnapPy for development purposes.  Just do

python dev/ccache.py command [options]

as you would with "setup.py".  For now, assumes Unix-like system with
gcc (not clang) as the compiler and relies on having "ccache" installed.

Unfortunately, no dependencies are deduced or checked here.  Do
"ccache -C" to clear the cache.  In particular, if you edit something
in "kernel" the corresponding "quad_double" code will not be
recompiled unless you manually touch said file.

---------

I tried to use the parallel compilation code given at:

http://stackoverflow.com/questions/11013851/speeding-up-build-process-with-distutils

without any success.
"""
import os, sys

# Make use ccache 
os.environ['CC'] = 'ccache gcc'
os.environ['CXX'] = 'ccache gcc++'

# Run the usual setup
sys.path = [os.path.abspath(os.curdir)] + sys.path
import setup
