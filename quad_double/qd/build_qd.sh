BASE_DIR=`pwd`
VERSION=2.3
LONG_VERSION=2.3.23
LIBZ=libz.${LONG_VERSION}.dylib
SRC_DIR=qd-${LONG_VERSION}
SRC_ARCHIVE=qd-${LONG_VERSION}.tar.gz
URL=https://www.davidhbailey.com/dhbsoftware/${SRC_ARCHIVE}
HASH=96f5521f14e00b0f9dbc7c95f4250420
if ! [ -e ${SRC_ARCHIVE} ]; then
    curl -O ${URL}
fi
ACTUAL_HASH=`md5 -q ${SRC_ARCHIVE}`
if [[ ${ACTUAL_HASH} != ${HASH} ]]; then
    echo Invalid hash value for ${SRC_ARCHIVE}
    exit 1
fi
if ! [ -d ${SRC_ARCHIVE} ]; then
    tar xvfz ${SRC_ARCHIVE}
fi
cd ${SRC_DIR}
mkdir -p lib
rm -rf lib/*
export CFLAGS="-arch x86_64 -mfpmath=sse"
export CXXFLAGS="-arch x86_64 -mfpmath=sse"
./configure --prefix ${BASE_DIR}/lib/intel --host x86_64-apple-darwin
make install
make clean
rm -rf config.log config.status
export CFLAGS="-arch arm64"
export CXXFLAGS="-arch arm64"
./configure --prefix ${BASE_DIR}/lib/arm --host arm64-apple-darwin
make install
make clean
rm -rf config.log config.status
cd ..
lipo -create lib/arm/lib/libqd.a lib/intel/lib/libqd.a -output lib/libqd.a
lipo -create lib/arm/lib/libqd.dylib lib/intel/lib/libqd.dylib -output lib/libqd.dylib

