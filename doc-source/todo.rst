To Do List
==========

- GUI
  
  - Allow you to get the info (position, radius) of a horoball by clicking it.  

- Documentation

  - More hyperlinks

  - Expand tutorial 

  - Make more screencasts

- snappy 
  
  - Longitude detection in manifolds with a single torus cusp (ie, the slope that dies in homology).  Requested by Saul.  

- Kernel 

  - Merge Jeff's changes into SnapPy's copy of the kernel.

  - After this, consider making SnapPy's copy the canonical one by
    definition. 
 
  - Improve the situation with uFatalErrors.  

- Minor 

  - dual_curves should really cache it's result and have this used by
    drill
  
  - One should be able to convert a SymmetryGroup to a Sage permutation group.   

  - Also, the SymmetryGroup presentation function should be wrapped.
    There is code for this in the old SnapPeaPython.  

- Ambitious

  - A new basic (sub)class: S3Knot (and/or S3Link).
 
  - Consider adding a HeegaardSplitting class 

  - Consider merging our t3m project and normal surface code into
    SnapPy. 

  - Redo much of Snap in the context of Sage/SnapPy.   

     - Add a method for computing tetrahedron shapes to arbitrary precision.

     - Add methods for computating invariant trace fields and related number
       fields.

     - Add a method which implements and extends Harriet Moser's
       algorithm to allow SnapPy to prove that a manifold is hyperbolic.


Development Basics
================================================

Submitting patches
-----------------------------------------


We're using Mercurial, and you can get a copy of the repository of
SnapPy and it's various components via::

   hg clone http://t3m.computop.org/hg/SnapPy
   hg clone http://t3m.computop.org/hg/plink
   hg clone http://t3m.computop.org/hg/Spherogram
   hg clone http://t3m.computop.org/hg/CyPari

After editing the files, commit your changes to the local repository via::

   hg commit -m "Fixed cache issue with sending from plink"

Then do::

   hg log -l 5
  
   changeset: 396:a5a0809a371d
   tag: tip
   user: Nathan Dunfield <nathan@dunfield.info>
   date: Fri Oct 02 	09:07:33 2009 -0500
   summary: Fixed cache issue with sending from plink
   ...

to get the changeset number(s) of your commits and then do::

  hg export -g 396:a5a0809a371d > plink_cache.patch
	
and mail us the file "plink_cache.patch".  


OS X
---------------------------

Here is how to get a clean development setup under OS X, versions
10.5-10.8.  

- Install Active Tcl/Tk 8.6 from `ActiveState
  <http://www.activestate.com/activetcl/>`_.

- Install Python 2.7.5 from Python.org using the `Mac Installer Disk Image
  <http://www.python.org/download/>`_.  There are currently two
  versions, one for 10.3 and up (ppc/i386) and one for 10.6 and up
  (i386/x86_64).  Either should work ok, but probably you want the
  second one.  Set your path so that "python" is::
      
    /Library/Frameworks/Python.framework/Versions/2.7/bin/python

- Then install `setuptools
  <https://pypi.python.org/pypi/setuptools/>`_, a Python
  package manager::

    curl -O http://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py
    python ez_setup.py 

  and use it to install the following packages::

    python -m easy_install mercurial   # Source code control software
    python -m easy_install Cython      # Used for Python-C interfacing
    python -m easy_install Sphinx      # For building the documentation
    python -m easy_install ipython     # Improved Python shell
    python -m easy_install py2app      # For making app bundles

- Get the source code from the repository, using the version of "hg" that
  is in the same directory as Python 2.7::

    hg clone http://t3m.computop.org/hg/plink
    hg clone http://t3m.computop.org/hg/Spherogram
    hg clone http://t3m.computop.org/hg/SnapPy

- Test the stand-alone link editor::

    cd plink
    python setup.py install
    python -m plink.app   # Link editor appears!

  This last command runs the script "plink/app.py"; the real code for
  the link editor is in "plink/__init__.py".

  To make sure it's using the right Tk, select "File->About Python..."
  and make sure the version is 8.6, not 8.4. or 8.5.  If it's an older
  version, go into "SnapPy/release_tools/tkinter-versions" and run the script
  "./install_tkinter 8.6".  (If you don't have both Python 3.2
  and 2.7 installed on your system, it will complain. But you can ignore
  this.)

  Building the proper Mac application bundle (not necessary for
  testing, so you can just skip this)::

    cd plink-app
    python setup.py py2app 
    open dist    # This directory contains the clickable app.  

- Now build SnapPy itself.  One builds it twice to generate the
  documentation, much of which is extracted from the installed module::

    cd ../SnapPy
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
    
Some major parts of the SnapPy codebase are:

- "SnapPy.pyx": The Cython interface to the SnapPea kernel
- "opengl/CyOpenGL*.pyx": The Cython interface to OpenGL*
- "snappy/app.py": The core GUI code
- "snappy/polyviewer.py": The GUI code for Dirichlet domains
- "snappy/horoviewer.py": The GUI code for horoball pictures
- "snappy/database.py": Interacts with the sqlite3 manifold database

In addition, Jeff's old prototype for a Tk-based UI can be found in
"misc/JeffsOldUI/SnapPeaGUI.py"; just run Python on this file to try it
out, after installing `PythonMegaWidgets <http://pmw.sf.net>`_.

Windows XP
-------------------------------------------------

Install `Python 2.7 <http://python.org>`_, `MinGW (including
g++, MSYS-base, and the MinGW Development Tookit) <http://mingw.org/wiki/Getting_Started>`_,
`Inno Setup <http://jrsoftware.org>`_, `Mercurial
<http://mercurial.selenic.com/downloads/>`_, and `PyReadine
<https://launchpad.net/pyreadline/+download>`_ via their binary
installers.  Due to `this bug  <http://bugs.python.org/issue12641>`_,
you need to edit by hand the file::

    c:Python27/Lib/distutils/cygwinccompiler.py

Inside the Mingw32CCompiler class there's a call to
"self.set_executables" and there you should remove all of the
"-mno-cygwin" options.  

Then install setuptools just by downloading `ez_setup.py
<http://peak.telecommunitycom/dist/ez_setup.py>`_ and double-clicking
it.  Then download the latest version of `Cython <http://cython.org>`_
into the directory "c:Python27".  In MSYS do the following::

   cd c:Python27
   tar xfz Cython-*.tar.gz
   cd Cython-*
   ../python.exe setup.py build -c mingw32
   ../python.exe setup.py install
   cd ../
   python.exe -m easy_install sphinx
   hg clone static-http://www.math.uic.edu/t3m/hg/SnapPy
   cd SnapPy
   sh build_pari.sh
   ../python.exe setup.py build -c mingw32
   ../python.exe setup.py install 
   ../python.exe setup.py build_docs
   ../python.exe setup.py install 
   cd ../
   python.exe -m snappy.app

If that works, install `py2exe <http://www.py2exe.org/>`_ via the binary installer.  Then::
 
   cd SnapPy/SnapPyExe
   export PATH=$PATH:/c/Python27:/c/Program\ Files/Inno\ Setup\ 5/

Now replace line 13 of make.py with the commented-out line 12.  Then::

  python make.py 
   

   
   






   
    
   
