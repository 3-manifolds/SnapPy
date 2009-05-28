#! /usr/bin/env python
#
# Script for building Togl under OS X.
#
# IMPORTANT NOTE: THIS DOESN'T SEEM TO ACTUALLY WORK.  

import os, sys, re


# /System/Library/FrameworksTk:
#TkVer = '8.4'
#SDK= '/Developer/SDKs/MacOSX10.5.sdk/System/Library'

# /Library/Frameworks Tk:
TkVer = '8.4'
SDK= '/Developer/SDKs/MacOSX10.5.sdk/Library'


tk_include_path = "-I" + SDK + "/Frameworks/Tk.framework/Versions/" + TkVer + "/Headers"
tk_includes = tk_include_path + " " + tk_include_path + "/tk-private"
tcl_include_path = "-I" + SDK + "/Frameworks/Tcl.framework/Versions/" + TkVer + "/Headers"
tcl_includes = tcl_include_path + " " + tcl_include_path + "/tcl-private"

# Can't use urllib or curl because of how sf.net works
if not os.path.exists('Togl2.0-src.tar.gz'):
    os.system('wget -nd http://downloads.sourceforge.net/togl/Togl2.0-src.tar.gz')

# Make sure we start with a clean copy.  
if os.path.exists('Togl2.0'):
    os.system('rm -rf Togl2.0')

os.system('tar xfz Togl2.0-src.tar.gz')
os.chdir('Togl2.0')
snappy_install="`pwd`/../snappy/darwin-tk" + TkVer
os.system('./configure --prefix=%s  --libdir=%s' % (snappy_install, snappy_install))
command ="make CFLAGS='-arch ppc -arch i386' CPPFLAGS='-arch ppc -arch i386 "  +   tk_includes + " " + tcl_includes + "'"
os.system(command)
os.system('make install-lib-binaries')
