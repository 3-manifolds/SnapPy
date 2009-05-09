#! /bin/bash
# This builds a fat pari library for OS X > 10.3.
# The includes that are installed are those for i386,
# which are definitely different from those for ppc
# (especially the assembly language in pariinl.h),
# but we don't use the parts where they differ.
#
# Eventually we need a switch for choosing the host.
# Probably we can get some help from distutils, but
# I need to understand eggs a little better to see
# how to do that.
curl --remote-name http://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.3.4.tar.gz
tar xvzf pari-2.3.4.tar.gz
cd pari-2.3.4
./Configure --prefix=`pwd` --host=ppc-darwin
cd Odarwin-ppc
make CFLAGS='-arch ppc' install-lib-sta
mv ../lib/libpari.a ../lib/ppc-libpari.a
cd ..
./Configure --prefix=`pwd` --host=i386-darwin
cd Odarwin-i386
make CFLAGS='-arch i386' install-lib-sta
mv ../lib/libpari.a ../lib/i386-libpari.a
lipo ../lib/ppc-libpari.a ../lib/i386-libpari.a -create -output ../lib/libpari.a
make install-include
ranlib ../*.a

