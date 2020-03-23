#! /usr/bin/env python
#
# Easy testing of the full Windows app without messing with the 
# installer.

import os
os.system('rm -rf build dist')
os.system('pyinstaller SnapPy_dbg.spec')
os.system('start dist/SnapPy_debug/SnapPy.exe')
