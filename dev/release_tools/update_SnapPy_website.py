#! /usr/bin/env python
from __future__ import print_function
import os, sys, re, glob, pkg_resources

web = '/afs/math.uic.edu/www/t3m/SnapPy/'
nest = '/afs/math.uic.edu/www/t3m/SnapPy-nest/'
nmd = '/afs/math.uic.edu/user/nmd'
hg_repo = '/afs/math.uic.edu/www/t3m/hg/SnapPy'

sandbox = nmd + "/SnapPy/"
frameworks = nmd +'/frameworks/'

print("Getting latest source code...")
os.chdir(sandbox)
os.system("hg pull") 
os.system("hg up")
print("Creating tarball...")
os.system("hg archive -t tar %s/SnapPy.tar" % nest)
os.system("gzip -f %s/SnapPy.tar" % nest)
os.chdir(nmd)

print("Deleting current egg...")
try:
    os.system("rm -Rf " + pkg_resources.get_distribution('snappy').location)
except pkg_resources.DistributionNotFound:
    pass 
    
print("Installing latest egg...")
os.system("python -m easy_install -U -f " + nest+ " snappy ")
pkg_resources.get_distribution('snappy').activate()
import snappy
version = snappy.version()
print("Current version appears to be: " + version)
docs = os.path.dirname(snappy.__file__) + os.sep + "doc"
print(docs)

print("Copying the documentation to the web...")
os.system("rm -rf " + web) 
os.system("cp -R " + docs + " " + web) 
os.system("ln -s " + nest + " " + web + "/get") 
os.system("ln -s " + web + " " + web + "/doc")

print("Here's all the matching eggs")
p = os.popen("ls -lat " + nest + "snappy-" + version + "*" + " " + nest + "InstallSnapPy.exe " + nest + "SnapPy.dmg " + nest + "SnapPy.tar.gz")
print(p.read())

ans = raw_input("Update 'current.txt' to " + version + " [N/y]: ")
if ans[:1] == 'y':
    open(nest + "current.txt", "w").write(version + "\n")
print("Contents of 'current.txt': " + open(nest+"current.txt").read().strip())

