"""
Build a temporary virtual environment and install SnapPy into it.
"""

import os, sys, re, urllib 

# If we don't have the virtual env script, go get it.

if not os.path.exists('virtualenv.py'):
    url = 'https://raw.github.com/pypa/virtualenv/master/virtualenv.py'
    try:
        import urllib
        urllib.urlretrieve(url, 'virtualenv.py')
    except AttributeError:
        import urllib.request
        urllib.request.urlretrieve(url, 'virtualenv.py')

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
        os.system(mypy + ' -m easy_install -f http://www.math.uic.edu/t3m/temp  ' + mod_name)
    os.system(mypy + ' ' + test_code) 

os.chdir(vir_env_dir)
test_module('SnapPy', '-m snappy.test')
test_module('CyPari', '-m cypari.test', False)
test_module('Spherogram', '-m spherogram.links.test', False)
