#!/bin/bash
UNAME=`uname`
if [ ! -e "qd-2.3.14.tar.gz" ] ; then
    curl -O -L http://crd-legacy.lbl.gov/~dhbailey/mpdist/qd-2.3.14.tar.gz
fi
if [ ! -e "qd-2.3.14" ]
then
    tar xvfz qd-2.3.14.tar.gz
fi
cd qd-2.3.14/include/qd
patch < ../../../qd.patch
cd ../..
if [ $UNAME == "Darwin" ] ; then
  ./configure --prefix=`pwd`/../../qd FCFLAGS='-m64' CFLAGS='-O3 -arch x86_64 -msse2 -mfpmath=sse' CXXFLAGS='-O3 -arch x86_64 -msse2 -mfpmath=sse'
  make install
elif [ $UNAME == "Linux" ] ; then
  ./configure --prefix=`pwd`/../../qd FCFLAGS='-m64' CFLAGS='-O3 -fPIC -msse2 -mfpmath=sse' CXXFLAGS='-O3 -fPIC -msse2 -mfpmath=sse'
  make install
else
    echo "Only tested on OS X and Linux so far."
fi