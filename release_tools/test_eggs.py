"""
Build a temporary virtual environment and install SnapPy into it.
"""

import os, sys, re

# Create virtual environment.

if sys.platform.startswith('win'):
    vir_env_dir = 'c:\SnapPyTestVE'
else: 
    vir_env_dir = '/tmp/SnapPyTestVE'
os.system('rm -rf ' + vir_env_dir)
#Note:  --no-site-packages  is implicit here.
os.system(sys.executable + ' virtualenv.py --distribute ' + vir_env_dir)
bindir = 'bin' if not sys.platform.startswith('win') else 'Scripts'
mypy = vir_env_dir + os.sep + bindir + os.sep + 'python'
def test_module(mod_name, test_code, install=True):
    print( '\n\n' + 30*'*' + '\n* ' + mod_name + '\n' + 30*'*' + '\n')
    if install:
        os.system(mypy + ' -m easy_install -f http://snappy.computop.org/get,http://dunfield.info/temp ' + mod_name)
    os.system(mypy + ' ' + test_code) 

os.chdir(vir_env_dir)
test_module('SnapPy', '-m snappy.test')
test_module('CyPari', '-m cypari.test')
test_module('Spherogram', '-m spherogram.links.test')
