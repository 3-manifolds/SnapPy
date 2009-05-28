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
./Configure --prefix=`pwd` --host=i386-linux
cd Olinux-i386
make install-lib-sta
make install-include

