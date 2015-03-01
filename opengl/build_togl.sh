#!/bin/bash
# SnapPy requires Togl, a OpenGL widget for the Tk libary.  For OS X
# and Linux, binaries for Togl are provided in the snappy/ directory,
# this script is only needed if, for some reason, those fail to work.
#
# This script builds this library in these cases:
# 
#   -Mac OS 10.5, against Tk 8.4
#   -Mac OS >= 10.5, against Tk 8.5.11 and Tk 8.6, using Togl2.1
#   -Linux, against Tk 8.5
# 

# You should be able to modify it, or follow the steps laid out here
# by hand, to compile Togl.  You will need to have the header files
# for Tcl and Tk available; for instance on Debian or Ubuntu you want
# the packages tcl-dev and tk-dev.  On some systems (e.g. Ubuntu 8.10),
# you may also need to install Xmu-dev and related packages.
#
#
# First, download the source if we don't have it already. Because
# Sourceforge is a little odd on how it does this, we use "wget"
# instead of "curl" or "lynx".  You may find it easier just to
# download the needed file with your webrowser..

set -e

if [ ! -e Togl2.0-src.tar.gz ]; then  
    echo "Downloading Togl2.0..."
    wget -nd http://downloads.sourceforge.net/togl/Togl2.0-src.tar.gz
fi
echo "Untaring Togl..."
tar xfz Togl2.0-src.tar.gz 
cd Togl2.0

# Set where we want Togl to be installed; this should be changed
# depending on your OS and version of Tk.  


if [ "$(uname)" = "Darwin" ] ; then  # If this is Mac OS X
    export SNAPPY_INSTALL=`pwd`/../../snappy/togl/darwin-tk8.4
else # Assume it's Linux
#   export SNAPPY_INSTALL=`pwd`/../snappy/linux2-tk8.5
   export SNAPPY_INSTALL=`pwd`/../../snappy/togl/linux2-tk8.6
fi

# Now build Togl. To configure Togl, the key is to find where
# tclConfig.sh and tkConfig.sh live.


if [ "$(uname)" = "Darwin" ] ; then  # If this is Mac OS X

./configure --prefix=$SNAPPY_INSTALL --libdir=$SNAPPY_INSTALL
make CFLAGS='-arch ppc -arch i386' \
CPPFLAGS='-arch ppc -arch i386 \
-I/Developer/SDKs/MacOSX10.5.sdk/System/Library/Frameworks/Tk.framework/Versions/8.4/Headers \
-I/Developer/SDKs/MacOSX10.5.sdk/System/Library/Frameworks/Tk.framework/Versions/8.4/Headers/tk-private \
-I /Developer/SDKs/MacOSX10.5.sdk/System/Library/Frameworks/Tcl.framework/Versions/8.4/Headers \
-I /Developer/SDKs/MacOSX10.5.sdk/System/Library/Frameworks/Tcl.framework/Versions/8.4/Headers/tcl-private'

else #Presume we have Linux here:
# uncomment to suit your system
#./configure --prefix=$SNAPPY_INSTALL --libdir=$SNAPPY_INSTALL \
#   --with-tcl=/usr/share/tcltk/tcl8.5/ --with-tk=/usr/share/tclk/tk8.5/ 
./configure --prefix=$SNAPPY_INSTALL --libdir=$SNAPPY_INSTALL \
   --with-tcl=/usr/lib/ --with-tk=/usr/lib/ 
#./configure --prefix=$SNAPPY_INSTALL --libdir=$SNAPPY_INSTALL \
#   --with-tcl=/usr/lib64/ --with-tk=/usr/lib64/ 
make

fi 


# and install it:

make install-lib-binaries

# Now we build Togl2.1, using the source code in Togl2.1-SnapPy.tgz

tar xvfz Togl2.1-SnapPy.tgz
cd Togl2.1
make -f Makefile.SnapPy
mv darwin-tk8.5 darwin-tk8.6 ../../snappy/togl 
