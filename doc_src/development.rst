Development Basics
================================================

Submitting patches
-----------------------------------------

The source code for SnapPy and its various parts are `hosted on
bitbucket <https://bitbucket.org/t3m>`_ as `Mercurial repositories
<http://mercurial.selenic.com>`_.   To contribute a patch, create a
free bitbucket account, fork the appropriate repository, and then send
us a pull request, as described in `this tutorial <https://confluence.atlassian.com/display/BITBUCKET/Fork+a+Repo%2C+Compare+Code%2C+and+Create+a+Pull+Request>`_.


OS X
---------------------------

Here is how to get a clean development setup under OS X.

- Install Active Tcl/Tk 8.6 from `ActiveState
  <http://www.activestate.com/activetcl/>`_.

- Install the latest Python 2.7 from Python.org using the `Mac
  Installer Disk Image <http://www.python.org/download/>`_.  There are
  currently two versions, one for 10.3 and up (ppc/i386) and one for
  10.6 and up (i386/x86_64); you want the second one.  Set your path
  so that "python" is::
      
    /Library/Frameworks/Python.framework/Versions/2.7/bin/python

- Python 2.7.9 and newer include `pip
  <https://pip.pypa.io/en/latest/index.html>`_ so use it upgrade and
  install the following packages::

    python -m pip install --upgrade setuptools
    python -m pip install virtualenv
    python -m pip install Cython      # Used for Python-C interfacing
    python -m pip install Sphinx      # For building the documentation
    python -m pip install ipython     # Improved Python shell
    python -m pip install py2app      # For making app bundles
    python -m pip install mercurial   # Source code control software 

- Get the source code from the repository.  The program "hg" was
  installed in the last step and lives in the same directory as Python 2.7::

    hg clone https://bitbucket.org/t3m/PLink
    hg clone https://bitbucket.org/t3m/Spherogram
    hg clone https://bitbucket.org/t3m/CyPari
    hg clone https://bitbucket.org/t3m/SnapPy

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

These instructions have been tested on Windows 7, 8 and 10 and quite
possibly work on XP and Vista as well. 

- Install `Python 2.7 <http://python.org>`_, specifically the 32 bit 
  version (Windows x86 not Windows x86-64) and The below instructions
  were checked with Python 2.7.11 and Inno Setup 5.5.4.  

- Install `MinGW (including g++, MSYS-base, and the MinGW Development
  Toolkit) <http://mingw.org/wiki/Getting_Started>`_, and open an MSYS
  terminal shell, which is where all the rest of the work will take
  place.  Alternatively, you can build everything except CyPari with
  this `Python-specific free version of Microsoft Visual C++
  <http://www.microsoft.com/en-us/download/details.aspx?id=44266>`_
  compiler.  If you would like to make your own installer you will
  also need `Inno Setup <http://www.jrsoftware.org/isdl.php>`_.

- If you wish to use the MinGW compiler to build everything, create
  a file "/c/Python27/Lib/distutils/distutils.cfg" consisting of::

    [build]	
    compiler=mingw32

  This tells Python to use the MinGW compilers to build all packages.
  You should skip this step if you're using the MSVC compiler instead.

- Make it so that MinGW, Python, and Inno Setup are all in
  your PATH by adding the below lines to the file "~/.profile"::

    PATH=/c/Python27:/c/Python27/Scripts:/c/mingw/bin:$PATH
    PATH=$PATH:'/c/Program\ Files\ \(x86\)/Inno\ Setup\ 5
    export PATH

- Python 2.7.9 and newer include `pip
  <https://pip.pypa.io/en/latest/index.html>`_ so let's use it
  to install the needed packages.::
  
    pip install -U pip      # Upgrades pip to the current version.
    pip install pyreadline 
    pip install sphinx
    pip install cython
    pip install ipython
    pip install pyx==0.12.1
    pip install mercurial   # Installs "hg", used in next step

- Fetch the latest development versions of the source straight from
  the t3m repository::

    hg clone https://bitbucket.org/t3m/plink
    hg clone https://bitbucket.org/t3m/Spherogram
    hg clone https://bitbucket.org/t3m/CyPari
    hg clone https://bitbucket.org/t3m/SnapPy

- Build and install each component of SnapPy, starting with CyPari::

    cd CyPari
    sh build_pari.sh
    python setup.py install

  If you elected to use the MSVC compiler to build SnapPy you must still
  use mingw32 to build CyPari; in this case the last command should be
  replaced by::

    python setup.py build --compiler=mingw32
    python setup.py

  Next build the other components::
  
    cd ../Spherogram
    python setup.py install
    cd ../plink 
    python setup.py install
    cd ../SnapPy
    python setup.py install
    cd ../

  Finally, start up the SnapPy app::
  
    python -m snappy.app 

- If that works, install `py2exe <http://www.py2exe.org/>`_ via the binary installer.  Then::

    cd SnapPy/windows_exe
    python make.py 

  builds the binary installer "InstallSnapPy.exe" for SnapPy.  

