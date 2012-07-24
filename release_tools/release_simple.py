"""
General script for building binary eggs for a t3m Python module.
"""

import os, sys, re, glob, subprocess, time

def run(command):
    os.system(command)


# Update (or get) the module source code from the repository
 
module_path = sys.argv[1]
if not os.path.exists(module_path):
    print('Getting module source...')
    run('hg clone http://math.uic.edu/t3m/hg/' + module_path)
os.chdir(module_path)
print('Updating module source...')
run('hg pull')
run('hg up')

# Build and install the module into the Python that called this script

if sys.platform.startswith('win'):
    run( sys.executable + ' setup.py build -c mingw32')
run( sys.executable + ' setup.py install' )

# For OS X, we need fake PPC/intel "fat" eggs.  

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
input('Hit any key when ready to begin copying to t3m:')
os.system('scp -p ' + ' '.join(eggs) + ' nmd@shell.math.uic.edu:t3m_web/SnapPy-nest')





