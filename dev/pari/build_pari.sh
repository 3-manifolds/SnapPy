#! /bin/bash
#
# This builds a fat pari library for OS X > 10.3 or a normal binary
# for Linux or Windows.  In the case of OS X, the includes that are 
# installed are those for i386, which are definitely different from
# those for ppc (especially the assembly language in pariinl.h), but
# we don't use the parts where they differ.  Eventually we need a
# switch for choosing the host.  Probably we can get some help from
# distutils, but I need to understand eggs a little better to see how
# to do that.
# For Windows (MinGW) some patches are applied.  These work around
# some UNIX functions that don't exist in MinGW (e.g. getuid()) and
# also work around the missing symbols (esp. time and localtime) that
# are listed in libmsvcr90.a but are not exported from msvcr90.dll .

if [ ! -e pari-2.5.1.tar.gz ]; then
    echo "Downloading Pari 2.5.1..."
    curl --remote-name http://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.5.1.tar.gz
fi
echo "Untarring Pari..."
tar xzf pari-2.5.1.tar.gz
cd pari-2.5.1

echo "Building Pari libary..." 
if [ "$(uname)" = "Darwin" ] ; then  # OS X
# 
#   export CFLAGS='-arch ppc -mmacosx-version-min=10.4 '
#    ./Configure --prefix=`pwd` --host=ppc-darwin
#    cd Odarwin-ppc
#    make install-lib-sta
#    cd ..
#    mv lib/libpari.a lib/ppc-libpari.a
#
    export CFLAGS='-arch i386 -mmacosx-version-min=10.4 '
    ./Configure --prefix=`pwd` --host=i386-darwin
    cd Odarwin-i386
    make install-lib-sta
    make install-include
    cd ..
    mv lib/libpari.a lib/i386-libpari.a
    mv include include32
#
    export CFLAGS='-arch x86_64'
    ./Configure --prefix=`pwd` --host=x86_64-darwin
    cd Odarwin-x86_64
    make install-lib-sta
    make install-include
    cd ..
    mv lib/libpari.a lib/x86_64-libpari.a
    mv include include64
#
#    lipo lib/ppc-libpari.a lib/i386-libpari.a lib/x86_64-libpari.a -create -output lib/libpari.a
    lipo lib/i386-libpari.a lib/x86_64-libpari.a -create -output lib/libpari.a
    ranlib lib/*.a
    ln -s include32 include
    #echo "******************* NOTICE! *********************"
    #echo "Pari is configured to build CyPari for 32-bit architectures."
    #echo "Change the symlink pari-2.3.4/include if you want to build for x86_64."
elif [ "$(uname)" = "MINGW32_NT-6.1" ] ; then # MinGW on Windows
    patch -p1 < ../mingw-pari.patch
    ./Configure --prefix=`pwd` --host=i386-mingw
    cd Omingw-i386
    make install-lib-sta
    make install-include
else  # Linux
    ./Configure --prefix=`pwd` 
    cd Olinux-*
    make install-lib-sta
    make install-include
fi 
