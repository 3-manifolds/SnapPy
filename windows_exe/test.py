#! /usr/bin/env python
#
# Easy testing of the full Windows app without messing with the 
# installer.

import os
os.system('rm -rf build dist')
os.system('pyinstaller SnapPy_py2.spec')
os.system('start dist/SnapPy/SnapPy.exe')
