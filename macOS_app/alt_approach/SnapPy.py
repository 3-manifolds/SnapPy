import os
import IPython
from snappy.app import main

# make sure we have a home directory where IPython can scribble
if os.name == 'nt' and 'HOME' not in os.environ:
    os.environ['HOME'] = os.environ['USERPROFILE']

def get_home_dir():
    return os.environ['HOME']

def get_ipython_dir():
    return os.path.join(os.environ['HOME'], '.ipython')

IPython.utils.path.get_home_dir = get_home_dir
IPython.utils.path.get_ipython_dir = get_ipython_dir()

main()
