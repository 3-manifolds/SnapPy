"""
Build a temporary virtual environment and install SnapPy into it.
"""

import os, sys, re

# Create virtual environment. 
vir_env_dir = '/tmp/SnapPyTestVE'
os.system('rm -rf ' + vir_env_dir)
os.system(sys.executable + ' virtualenv.py --distribute --no-site-packages ' + vir_env_dir)
mypy = vir_env_dir + '/bin/python'

def test_module(mod_name, test_code, install=True):
    print( '\n\n' + 30*'*' + '\n* ' + mod_name + '\n' + 30*'*' + '\n')
    if install:
        os.system(mypy + ' -m easy_install -f http://snappy.computop.org/get ' + mod_name)
    os.system(mypy + ' ' + test_code) 

os.chdir('/tmp')
test_module('SnapPy', '-m snappy.test')
test_module('CyPari', '-m cypari.test')
test_module('Spherogram', '-m spherogram.links.test')
