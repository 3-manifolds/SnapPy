#! /usr/bin/env python

import os, sys, re, glob

# Now put it on the webpage:

user = os.environ['USER']
if user in ['nmd', 'dunfield']:
    print "Hi there Nathan..."
    address = "nmd@shell.math.uic.edu"
if user == 'culler':
    print "Hi there Marc..."
    address = "culler@threlfall.math.uic.edu"

for file in glob.glob("../dist/*-intel.egg"):
    copy = file.replace("-intel", "-fat")
    os.system("cp " + file + " " + copy)
    
os.system("chmod g+w dmg-maker/SnapPy.dmg ../dist/*.egg")
os.system("scp -p dmg-maker/SnapPy.dmg ../dist/*.egg %s:/afs/math.uic.edu/www/t3m/SnapPy-nest" % address)
