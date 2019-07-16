We currently relases an app and wheels for 10.6 (or is it 10.9?) and up
for 64-bit intel only.



Building Tcl/Tk from source
---------------------------

This isn't typically necessary, but sometimes it's necessary to stay
on the bleeding edge.  Initial installation::

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

  cd /pkgs/tcl_tk

  cd tcl; fossil update core-8-6-branch; make -C macosx
  sudo make -C macosx install

  cd ../tk; fossil update core-8-6-branch; make -C macosx
  sudo make -C macosx install

Assuming one has installed a recent Python 3.7 (say) from python.org,
then you have to overide its internal copy of Tk as follows::

  sudo cp  /Library/Frameworks/Tk.framework/Versions/8.6/Tk \
           /Library/Frameworks/Python.framework/Versions/3.7/lib/libtk8.6.dylib
  sudo cp /Library/Frameworks/Tk.framework/Versions/8.6/libtkstub8.6.a \
          /Library/Frameworks/Python.framework/Versions/3.7/lib/
  sudo cp -R \
      /Library/Frameworks/Tk.framework/Versions/8.6/Resources/Scripts/ \
      /Library/Frameworks/Python.framework/Versions/3.7/lib/tk8.6/

  
