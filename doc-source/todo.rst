To Do List
==========

- GUI

  - Make the colors of link components in PLink match the colors
    used by SnapPy's horoball viewer. 

- snappy

- Documentation

  - More hyperlinks
  - Make screencast


- Tasks for another day:
  - Fix magic_ls in the interpreter.

  - Fix plink DT code issue for links and snap
   
  - Splittings (also requested by Saul)

  - Longitude detection in manifolds with a single torus cusp (ie, the slope that dies in homology).  Requested by Saul.  

  - dual_curves should really cache it's result and have this used by
    drill
  
  - Abelian group should he able to take any input and put it in
    canonical form, rather than simply insisting it be that way already. 
    (Cf  kernel_code/abelian_group.c/compress_abelian_group())

  - One should be able to convert a SymmetryGroup to a Sage permutation group.   

  - Also, the SymmetryGroup presentation function should be wrapped.
    There is code for this in the old SnapPeaPython.  

Development Basics
================================================

Submitting patches
-----------------------------------------


We're using mercurial, and you can get a copy of the repository via::

   hg clone static-http://www.math.uic.edu/t3m/hg/SnapPy

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
10.4, 10.5 and 10.6.  

- First, install Python 2.7 using the `Mac Installer Disk Image
  <http://http://www.python.org/download/>`_.  Be sure use the
  installer which works for OS X versions 10.3 and up, not the one for
  10.5 and up; SnapPy's graphical features may not work with the
  latter installer because of `this bug
  <http://bugs.python.org/issue9227>`_.   Set your path so that
  "python" is::
      
    /Library/Frameworks/Python.framework/Versions/2.7/bin/python

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

- Install Active Tcl/Tk 8.4.19 (not 8.5 or 8.6) from `ActiveState
  <http://www.activestate.com/activetcl/>`_.

- Get the source code from the repository, using the version "hg" that
  is in the same directory as Python 2.7::

    hg clone static-http://www.math.uic.edu/t3m/hg/plink
    hg clone static-http://www.math.uic.edu/t3m/hg/SnapPy

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

Windows XP
-------------------------------------------------

Install `Python 2.7 <http://python.org>`_, `MinGW (including
g++, MSYS-base, and the MinGW Development Tookit) <http://mingw.org/wiki/Getting_Started>`_,
`Inno Setup <http://jrsoftware.org>`_, `Mercurial
<http://mercurial.selenic.com/downloads/>`_, and `PyReadine
<https://launchpad.net/pyreadline/+download>`_ via their binary
installers.  Then install setuptools just by downloading `ez_setup.py
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
   

   
   






   
    
   
