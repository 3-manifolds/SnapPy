#!/bin/bash
#wget -nd http://downloads.sourceforge.net/togl/Togl2.0-src.tar.gz
tar xvfz Togl2.0-src.tar.gz 
cd Togl2.0
export SNAPPY_INSTALL=`pwd`/../snappy/linux2-tk8.5
./configure --prefix=$SNAPPY_INSTALL --libdir=$SNAPPY_INSTALL --with-tcl=/usr/share/tcltk/tcl8.5/ --with-tk=/usr/share/tcltk/tk8.5/ 
make
make install-lib-binaries