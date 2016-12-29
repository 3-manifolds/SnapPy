We currently relases an app and eggs for 10.5+ and intel (both 32 bit and 64 bit).  

The Python used to build things was compiled with:

./configure --enable-framework --with-framework-name=Python-10.5-intel --enable-universalsdk=/Developer/SDKs/MacOSX10.5.sdk --with-universal-archs=intel


Building Tcl/Tk from source
------------------------------

This isn't typically necessary, but sometimes it's necessary to stay on the bleeding edge.  Initial installation::

  brew install fossil    # Obscure DVCS
  mkdir /pkgs/tcl_tk
  cd /pkgs/tcl_tk
  fossil clone http://core.tcl.tk/tk tk.fossil
  fossil clone http://core.tcl.tk/tcl tcl.fossil
  mkdir tcl tk
  cd tcl
  fossil open ../tcl.fossil
  cd ../tk
  fossil open ../tk.fossil

Update and build::

  /pkgs/tcl_tk

  cd tcl; fossil update core-8-6-branch; make -C macosx
  sudo make -C macosx install

  cd ../tk; fossil update core-8-6-branch; make -C macosx
  sudo make -C macosx install

