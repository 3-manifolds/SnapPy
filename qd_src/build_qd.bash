#!/bin/bash
set -e
UNAME=`uname -s | cut -b -6`
QDVERSION="qd-2.3.14"
QDTARBALL="$QDVERSION.tar.gz"
QDURL="http://crd-legacy.lbl.gov/~dhbailey/mpdist/$QDTARBALL"

# grab the qd source, if necessary
if [ ! -e $QDTARBALL ] ; then
     echo "Downloading $QDVERSION..."
     python -c "import urllib; urllib.urlretrieve('$QDURL', '$QDTARBALL')"
fi

if [ ! -e  $QDVERSION ]
then
    echo "Untarring $QDTARBALL..."
    tar xfz $QDTARBALL
fi
# apply our patches
echo "Applying patches"
cd $QDVERSION/include/qd
patch < ../../../qd.patch
cd ../..
# build and install
echo "Building QD library"
if [ $UNAME == "Darwin" ] ; then
  ./configure --prefix=`pwd`/../../qd --enable-fortran=no FCFLAGS='-m64' CFLAGS='-O3 -arch x86_64 -msse2 -mfpmath=sse -mieee-fp' CXXFLAGS='-O3 -arch x86_64 -msse2 -mfpmath=sse -mieee-fp'
  make install
elif [ $UNAME == "Linux" ] ; then
  ./configure --prefix=`pwd`/../../qd --enable-fortran=no FCFLAGS='-m64' CFLAGS='-fPIC -msse2 -mfpmath=sse -mieee-fp' CXXFLAGS='-fPIC -msse2 -mfpmath=sse -mieee-fp'
  make install
elif [ $UNAME == "MINGW3" ] ; then
  ./configure --prefix=`pwd`/../../qd --enable-fortran=no CFLAGS='-O3 -msse2 -mfpmath=sse -mieee-fp -mstackrealign' CXXFLAGS='-O3 -msse2 -mfpmath=sse -mieee-fp -mstackrealign'
  cp /c/mingw/include/float.h .
  make install
elif [ $UNAME == "MINGW_" ] ; then
  ./configure --prefix=`pwd`/../../qd --enable-fortran=no CFLAGS='-O3 -msse2 -mfpmath=sse -mieee-fp -mstackrealign' CXXFLAGS='-O3 -msse2 -mfpmath=sse -mieee-fp -mstackrealign'
  cp /c/mingw/include/float.h .
  make install
else
    echo "Only tested on OS X, Linux, and MinGW so far."
fi
