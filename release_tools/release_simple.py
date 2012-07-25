"""
General script for building binary eggs for a t3m Python module.
"""

import os, sys, re, glob, subprocess, time

def run(command):
    os.system(command)


# Update (or get) the module source code from the repository
 
module_path = sys.argv[1]
if not os.path.exists(module_path):
    print('*** Getting module source...')
    run('hg clone http://math.uic.edu/t3m/hg/' + module_path)
os.chdir(module_path)
print('*** Updating module source...')
run('hg pull')
run('hg up')

# Build and install the module into all the standard Pythons

if sys.platform.startswith('linux'):
    pythons = ['/home/dunfield/python' + py + '/bin/python' for py in ['26', '27', '32']]
elif sys.platform.startswith('darwin'):
    pythons = ['/Library/Frameworks/Python.framework/Versions/2.6/bin/python',
               '/Library/Frameworks/Python-10.5-intel.framework/Versions/2.7/bin/python2.7',
               '/Library/Frameworks/Python-10.5-intel.framework/Versions/3.2/bin/python3.2']
elif sys.platform.startswith('win'):
    python = ['/c/Python26/python.exe', '/c/Python27/python.exe', '/c/Python32/python.exe']
               
for python in pythons:
    print "\n*** Building version with " + python
    if sys.platform.startswith('win'):
        run( python + ' setup.py build -c mingw32')
    run(python + ' setup.py install' )

# For OS X, we need fake PPC/intel "fat" eggs.  

if sys.platform.startswith('darwin'):
    for file in glob.glob("dist/*-intel.egg"):
        copy = file.replace("-intel", "-fat")
        os.system("cp " + file + " " + copy)

# Locate eggs that are less than a day old
    
cut_time = time.time() - 24*60*60
eggs = [egg for egg in glob.glob('dist/*.egg') if os.path.getmtime(egg) > cut_time]

if sys.platform.startswith('win'):
    eggs = [e.replace('\\', '/') for e in eggs]

try: input = raw_input
except: pass
input('***Hit any key when ready to begin copying to t3m***')
os.system('scp -p ' + ' '.join(eggs) + ' nmd@shell.math.uic.edu:t3m_web/SnapPy-nest')





