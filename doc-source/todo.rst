To Do List
==========

- GUI

- snappy

- Documentation

  - More hyperlinks
  - Make screencast

- Tasks for another day:
   
  - Splittings 

  - Starting with generators of a lattice in PSL(2,C), building a
    fundamental domain and getting a triangulation of the corresponding
    manifold.  (Useful for building arithmetic hyperbolic 3-manifolds.)

  - dual_curves should really cache it's result and have this used by
    drill
  
  - Abelian group should he able to take any input and put it in
    canonical form, rather than simply insisting it be that way already. 
    (Cf  kernel_code/abelian_group.c/compress_abelian_group())

  - One should be able to convert a SymmetryGroup to a Sage permutation group.   

  - Also, the SymmetryGroup presentation function should be wrapped.
    There is code for this in the old SnapPeaPython.  

  - Make it accept the alphanumeric DT codes that snap uses.


Development Basics
=================================

Here is how to get a clean development setup under OS X, either 10.4
or 10.5.   

- First, install Python 2.6 using the `Mac Installer Disk Image 
  <http://http://www.python.org/download/>`_.  Set your path so that
  "python" is::
      
    /Library/Frameworks/Python.framework/Versions/2.6/bin/python

- Then install `Setuptools
  <http://peak.telecommunity.com/DevCenter/setuptools>`_, a Python
  package manager::

    curl -O http://peak.telecommunity.com/dist/ez_setup.py
    python ez_setup.py  

  and use it to install the following packages::

    python -m easy_install mercurial   # Source code control software
    python -m easy_install Cython      # Used for Python-C interfacing
    python -m easy_install Sphinx      # For building the documentation
    python -m easy_install ipython     # Improved Python shell
    python -m easy_install py2app      # For making app bundles
    touch /Library/Frameworks/Python.framework/Versions/2.6/Resources/version.plist # Fixes bug in py2app

- Install Active Tcl/Tk 8.4.19 (not 8.5 or 8.6) from `ActiveState
  <http://www.activestate.com/activetcl/>`_.

- Get the source code from the repository, using the version "hg" that
  is in the same directory as Python 2.6::

    hg clone static-http://www.math.uic.edu/~t3m/hg/plink
    hg clone static-http://www.math.uic.edu/~t3m/hg/SnapPy

  Next, make it so that Python uses the version of Tk that was just
  installed, not the older version that ships with OS X::

    sudo cp SnapPy/tkinter-versions/_tkinter-8.4.19.so /Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/lib-dynload/_tkinter.so

- Test the stand-alone link editor::

    cd plink
    python setup.py install
    python -m plink.app   # Link editor appears!

  This last command runs the script "plink/app.py"; the real code for
  the link editor is in "plink/__init__.py".

  Building the proper Mac application bundle (not necessary for
  testing, typically)::

    cd plink-app
    python setup.py py2app 
    open dist    # This directory contains the clickable app.  

- Now build SnapPy itself.  One builds it twice to generate the
  documentation, much of which is extracted from the installed module::

    cd ../../SnapPy
    sh build_pari.sh     # Used to compute homology
    python setup.py install
    python setup.py build_docs install  

  If "." is in your path, you'll need to change directory before starting
  SnapPy; otherwise it will attempt to load "./snappy" which lacks the
  binary module::

    cd SnapPyApp
    python -m snappy.app   #SnapPy starts!

  To build the clickable app, just do the following in the SnapPyApp
  directory::

    python setup.py py2app
    
The some parts of the SnapPy codebase are:

- "SnapPy.pyx": The Cython interface to the SnapPea kernel
- "CyOpenGL*.pyx": The Cython interface to OpenGL*
- "snappy/app.py": The core GUI code
- "snappy/polyviewer.py": The GUI code for Dirichlet domains
- "snappy/horoviewer.py": The GUI code for horoball pictures

In addition, Jeff's old prototype for a Tk-based UI can be found in
"JeffsOldUI/SnapPeaGUI.py"; just run Python on this file to try it
out, after installing `PythonMegaWidgets <http://pmw.sf.net>`_.





   
    
   
