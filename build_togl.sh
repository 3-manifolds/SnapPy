#!/bin/bash
wget -nd http://downloads.sourceforge.net/togl/Togl2.0-src.tar.gz
tar xvfz Togl2.0-src.tar.gz 
cd Togl2.0
export SNAPPY_INSTALL=`pwd`/../SnapPy/darwin-tk8.4
./configure --prefix=$SNAPPY_INSTALL --libdir=$SNAPPY_INSTALL
make CPPFLAGS='-I/Developer/SDKs/MacOSX10.5.sdk/System/Library/Frameworks/Tk.framework/Versions/8.4/Headers \
-I/Developer/SDKs/MacOSX10.5.sdk/System/Library/Frameworks/Tk.framework/Versions/8.4/Headers/tk-private \
-I /Developer/SDKs/MacOSX10.5.sdk/System/Library/Frameworks/Tcl.framework/Versions/8.4/Headers \
-I /Developer/SDKs/MacOSX10.5.sdk/System/Library/Frameworks/Tcl.framework/Versions/8.4/Headers/tcl-private'
make install-lib-binaries