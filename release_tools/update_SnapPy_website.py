#! /usr/bin/env python
from __future__ import print_function
import os, sys, re, glob

web = '/afs/math.uic.edu/www/t3m/SnapPy/'
nest = '/afs/math.uic.edu/www/t3m/SnapPy-nest/'
nmd = '/afs/math.uic.edu/user/nmd'
hg_repo = '/afs/math.uic.edu/www/t3m/hg/SnapPy'

sandbox = nmd + "/SnapPy/"
frameworks = nmd +'/frameworks/'

print("Creating tarball...")
os.chdir(hg_repo)
os.system("hg archive -t tar %s/SnapPy.tar" % nest)
os.system("gzip -f %s/SnapPy.tar" % nest)

#print("Installing latest egg...")
#os.system("python -m easy_install -U -f http://snappy.computop.org/get  snappy ")
#import snappy.version
#version = snappy.version.version
#print("Current version appears to be: " + version)
#import snappy
#docs = os.path.dirname(snappy.__file__) + os.sep + "doc"

#print("Copying the documentation to the web...")
#os.system("rm -rf " + web) 
#os.system("cp -R " docs + " " + web) 
#os.system("ln -s " + nest + " " + web + "/get") 
#os.system("ln -s " + web + " " + web + "/doc")

#print("Here's all the matching eggs")
#os.popen("ls -lat " + nest + "/snappy-" + version + "*" + " " + nest + "/InstallSnapPy.exe " + nest + "/SnapPy.dmg " + nest + "SnapPy.tar.gz")

