#! /usr/bin/env python
#
# To just copy the result, add a '-c' argument.  

import os, sys, re

home = os.environ['HOME']
source_code = home + "/SnapPy/"
python26 = home + "/python26/bin/python"
python27 = home + "/python27/bin/python"
python32 = home + "/python32/bin/python"

if len(sys.argv) == 1:
    os.chdir(source_code)
    os.system("hg pull")
    os.system("hg up")
    os.system("rm dist/*.egg")
    for python in [python26, python27]: 
        os.system(python + " setup.py install")
        os.system(python + " setup.py build_docs install")
