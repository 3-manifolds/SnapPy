#!/bin/bash
UNAME=`uname -s | cut -b -6`
# grab the qd source, if necessary
if [ ! -e "qd-2.3.14.tar.gz" ] ; then
   if [ $UNAME == "Darwin" ] ; then
     curl -O -L http://crd-legacy.lbl.gov/~dhbailey/mpdist/qd-2.3.14.tar.gz
   else
     wget http://crd-legacy.lbl.gov/~dhbailey/mpdist/qd-2.3.14.tar.gz
   fi
fi
# unpack the archive, if necessary
if [ ! -e "qd-2.3.14" ]
then
    tar xvfz qd-2.3.14.tar.gz
fi
# apply our patches
cd qd-2.3.14/include/qd
patch < ../../../qd.patch
cd ../..
# build and install
if [ $UNAME == "Darwin" ] ; then
  ./configure --prefix=`pwd`/../../qd FCFLAGS='-m64' CFLAGS='-O3 -arch x86_64 -msse2 -mfpmath=sse' CXXFLAGS='-O3 -arch x86_64 -msse2 -mfpmath=sse'
  make install
elif [ $UNAME == "Linux" ] ; then
  ./configure --prefix=`pwd`/../../qd FCFLAGS='-m64' CFLAGS='-O3 -fPIC -msse2 -mfpmath=sse' CXXFLAGS='-O3 -fPIC -msse2 -mfpmath=sse'
  make install
elif [ $UNAME == "MINGW3" ] ; then
  ./configure --prefix=`pwd`/../../qd CFLAGS='-O3 -msse2 -mfpmath=sse -mstackrealign' CXXFLAGS='-O3 -msse2 -mfpmath=sse -mstackrealign'
  cp /c/mingw/include/float.h .
  make install
elif [ $UNAME == "MINGW_" ] ; then
  ./configure --prefix=`pwd`/../../qd CFLAGS='-O3 -msse2 -mfpmath=sse -mstackrealign' CXXFLAGS='-O3 -msse2 -mfpmath=sse -mstackrealign'
  cp /c/mingw/include/float.h .
  make install
else
    echo "Only tested on OS X, Linux and MinGW so far."
fi
