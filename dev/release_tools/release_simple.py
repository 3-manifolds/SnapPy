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
    run('hg clone https://bitbucket.org/t3m/' + module_path)
os.chdir(module_path)
print('*** Updating module source...')
run('hg pull')
run('hg up')

# Build and install the module into all the standard Pythons

if sys.platform.startswith('linux'):
    pythons = ['/home/dunfield/python' + py + '/bin/python' for py in ['26', '27']]
    commands = ['install']
elif sys.platform.startswith('darwin'):
    special_python = '/Library/Frameworks/Python-10.5-intel.framework/Versions/2.7/bin/python2.7'
    if os.path.exists(special_python):
        pythons = [special_python]
        commands = ['install', 'sdist']
    else:
        pythons = ['/Library/Frameworks/Python.framework/Versions/2.7/bin/python2.7']
        commands = ['bdist_wheel']
elif sys.platform.startswith('win'):
    pythons = ['c:\Python27\python.exe']
    commands = ['install', 'bdist_wheel']
               
for python in pythons:
    print "\n*** Building version with " + python
    if sys.platform.startswith('win'):
        run( python + ' setup.py build -c mingw32')
    for command in commands:
        run(python + ' setup.py ' + command)


# Locate eggs that are less than a day old
    
cut_time = time.time() - 24*60*60
eggs = [egg for egg in glob.glob('dist/*.egg') if os.path.getmtime(egg) > cut_time]

if sys.platform.startswith('win'):
    eggs = [e.replace('\\', '/') for e in eggs]

# For OS X, we need fake PPC/intel "fat" eggs.  
# 
#if sys.platform.startswith('darwin'):
#    for egg in eggs[:]:
#        copy = egg.replace("-intel", "-fat")
#        os.system("cp " + egg + " " + copy)
#        eggs.append(copy)


#try: input = raw_input
#except: pass
#input('***Hit any key when ready to begin copying to t3m***')
#os.system('scp -p ' + ' '.join(eggs) + ' nmd@shell.math.uic.edu:t3m_web/SnapPy-nest')





