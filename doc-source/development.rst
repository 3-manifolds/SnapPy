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

- Install the latest Python 2.7 from Python.org using the `Mac
  Installer Disk Image <http://www.python.org/download/>`_.  There are
  currently two versions, one for 10.3 and up (ppc/i386) and one for
  10.6 and up (i386/x86_64).  Either should work ok, but probably you
  want the second one.  Set your path so that "python" is::
      
    /Library/Frameworks/Python.framework/Versions/2.7/bin/python

- Then install `setuptools
  <https://pypi.python.org/pypi/setuptools/>`_, a Python
  package manager::

    curl -O https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py
    python ez_setup.py 

  and use it to install the following packages::

    python -m easy_install mercurial   # Source code control software
    python -m easy_install Cython      # Used for Python-C interfacing
    python -m easy_install Sphinx      # For building the documentation
    python -m easy_install ipython     # Improved Python shell
    python -m easy_install py2app      # For making app bundles

- Get the source code from the repository.  The program "hg" was
  installed in the last step and lives in the same directory as Python 2.7::

    hg clone http://t3m.computop.org/hg/plink
    hg clone http://t3m.computop.org/hg/Spherogram
    hg clone http://t3m.computop.org/hg/CyPari
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

- Build and install Spherogram::

    cd ../Spherogram
    python setup.py install

- Build and install CyPari::

    cd ../CyPari
    python setup.py install

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

- "SnapPy.pyx", "SnapPycore.pxi", "SnapPy.pxi": The Cython interface
  to the SnapPea kernel
- "opengl/CyOpenGL*.pyx": The Cython interface to OpenGL*
- "snappy/app.py": The core GUI code
- "snappy/polyviewer.py": The GUI code for Dirichlet domains
- "snappy/horoviewer.py": The GUI code for horoball pictures
- "snappy/database.py": Interacts with the sqlite3 manifold database

In addition, Jeff's old prototype for a Tk-based UI can be found in
"misc/JeffsOldUI/SnapPeaGUI.py"; just run Python on this file to try it
out, after installing `PythonMegaWidgets <http://pmw.sf.net>`_.

Windows
-------------------------------------------------

These instructions have been tested on Windows 7 and 8 and quite
possibly work on XP and Vista as well. 

- Install `Python 2.7 <http://python.org>`_, specifically the 32 bit 
  version (Windows x86 not Windows x86-64) and also `Inno Setup
  <http://jrsoftware.org>`_.  The below instructions were checked with
  Python 2.7.6 and Inno Setup 5.5.4.  

- Install `MinGW (including g++, MSYS-base, and the MinGW Development
  Toolkit) <http://mingw.org/wiki/Getting_Started>`_, and open an MSYS
  terminal shell, which is where all the rest of the work will take
  place. 

- Create a file "/c/Python27/Lib/distutils/distutils.cfg" consisting
  of::

    [build]	
    compiler=mingw32

  This tells Python to use the MinGW compilers.  

- Make it so that MinGW, Python, and Inno Setup are all in
  your PATH by adding the below lines to the file "~/.profile"::

    PATH=/c/Python27:/c/Python27/Scripts:/c/mingw/bin:$PATH
    PATH=$PATH:'/c/Program Files/Inno Setup 5'
    export PATH

- Install `"pip"
  <http://www.pip-installer.org/en/latest/installing.html>`_, which in
  turn installs both "setuptools" and "easy_install".  

- Install various Python packages::
  
	pip install pyreadline 
	pip install sphinx
	pip install cython
	pip install ipython
	pip install --allow-all-external pyx==0.12.1
	pip install mercurial   # Installs "hg", used in next step

- Fetch the latest development versions of the source straight from
  the repository::

        hg clone http://t3m.computop.org/hg/CyPari
	hg clone http://t3m.computop.org/hg/spherogram
	hg clone http://t3m.computop.org/hg/plink
	hg clone http://t3m.computop.org/hg/SnapPy

- Build and install each piece of the library in turn, and then start SnapPy::

    cd CyPari
    python setup.py install
    cd ../Spherogram
    python setup.py install
    cd ../plink 
    python setup.py install
    cd ../SnapPy
    python setup.py install
    cd ../
    python -m snappy.app 

- If that works, install `py2exe <http://www.py2exe.org/>`_ via the binary installer.  Then::

    cd SnapPy/SnapPyExe
    python make.py 

  builds the binary installer "InstallSnapPy.exe" for SnapPy.  

