Development Basics
================================================

Submitting patches
-----------------------------------------

The source code for SnapPy and its various parts are `hosted on
bitbucket <https://bitbucket.org/t3m>`_ as `Mercurial repositories
<https://www.mercurial-scm.org/>`_.  To contribute a patch, create a
free bitbucket account, fork the appropriate repository, and then send
us a pull request, as described in `this tutorial
<https://confluence.atlassian.com/bitbucket/create-a-pull-request-774243413.html>`_.


OS X
---------------------------

Here is how to get a clean development setup under OS X.

- Install the latest Python 3.7 from Python.org using the `Mac
  Installer Disk Image <http://www.python.org/download/>`_.  There are
  currently two versions, one for 10.6 and up (64-bit/32-bit) and one
  for 10.9 and up (64-bit); you want the second one.  Set your path so
  that "python3" is::
      
    /Library/Frameworks/Python.framework/Versions/3.7/bin/python3

- Use `pip <https://pip.pypa.io/en/latest/index.html>`_ to install the
  following packages::

    python3 -m pip install --upgrade setuptools virtualenv wheel pip
    python3 -m pip install cython        # Used for Python-C interfacing
    python3 -m pip install sphinx        # For building the documentation
    python3 -m pip install ipython       # Improved Python shell
    python3 -m pip install py2app        # For making app bundles
    python3 -m pip install pyx

- Get the source code from the repository, using Mercurial. For
  example you can install Mercurial via its `package installer
  <https://www.mercurial-scm.org/>`_.  Now do::

    hg clone https://bitbucket.org/t3m/plink
    hg clone https://bitbucket.org/t3m/spherogram
    hg clone https://bitbucket.org/t3m/snappy

- Test the stand-alone link editor::

    cd plink
    python3 setup.py pip_install
    python3 -m plink.app   # Link editor appears!

  This last command runs the script "plink/app.py"; the real code for
  the link editor is in "plink/__init__.py".

- Build and install Spherogram::

    cd ../spherogram
    python3 setup.py pip_install
    python3 setup.py test

- Now build SnapPy itself.  One builds it twice to generate the
  documentation, much of which is extracted from the installed module::

    cd ../snappy
    python3 setup.py pip_install
    python3 setup.py test   # Run the tests; pretty picture should appear.
    python3 setup.py build_docs pip_install
    python3 -m snappy.app   #SnapPy starts!

  To build the clickable app, just do the following::

    cd mac_osx_app
    python3 release.py --manual

  though for general development purposes ``python3 -m snappy.app`` is
  usually the way to go.
    
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
-------

These instructions have been tested on Windows 10, and describe
setting up a development environment using the (free) MSVC
compiler. To build the CyPari subcomponent, which few will want or
need to do, one must install additional tools.

- Install `Python 3.7 <https://www.python.org/downloads/windows/>`_,
  specifically the 32 bit version (Windows x86 not Windows x86-64).
  These instructions assume it has been installed in the directory
  ``C:\Python37`` which is not the default.

- With Python 3.7.4, you need the 2017 edition of the MSVC command
  line tools.  You can get them by using Microsoft's free `Build Tools for Visual Studio 2019
  <https://visualstudio.microsoft.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=16>`_
  installer and selecting the "C++ command line tools" and then in the
  sidebar selection only MCV 141. 

  
  If you would like to make your own installer you will also need
  `Inno Setup <http://www.jrsoftware.org/isdl.php>`_, specifically the
  unicode version; tested with version 5.5.9.

- Install whichever version of `MSYS2 <http://msys2.github.io>`_ is
  appropriate for your version Windows.  Most commonly, you will have
  a 64-bit Windows and hence want the "x86_64" installer; for
  concreteness the rest of these instructions assume this. (Technical
  note: even if you want to build 32-bit binaries, if your Windows is
  64-bit you want the x86_64 installer.) Follow the instructions on
  the webpage to update everything to the very latest MSYS2
  (``pacman -Sy pacman; pacman -Syu; pacman -Su`` etc.).

- Make a shortcut to ``c:\msys64\msys2.exe`` as you will be using it all
  the time; alternatively, pin ``mys2.exe`` to your taskbar.  

- Install some additional packages::

    pacman -S git make nano openssh perl tar unzip wget winpty patch

- Install your favorite text editor, for example you can install Emacs
  via::

    pacman -S  mingw-w64-x86_64-emacs

- Make it so that MinGW, Python, and Inno Setup are all in your PATH,
  as well as work around some stupid bug, by making the end of your
  "~/.bash_profile" file to read::

    PATH=/c/Python37:/c/Python37/Scripts:$PATH
    PATH=$PATH:'/c/Program Files (x86)/Inno Setup 5'
    export PATH
    alias emacs="/mingw64/bin/emacs"
    winpty bash; exit

  For example, do::

    nano ~/.bash_profile

- Make sure you have the right version of Python in your path by
  typing::

    python --version

  You should see something like ``Python 3.7.4``.

- Use pip to install some basic tools::
  
    python -m pip install --upgrade pip setuptools wheel  # Upgrades pip to the current version.
    python -m pip install pyreadline sphinx cython cypari pyx

- Install the source-code control system `Mercurial <https://www.mercurial-scm.org/>`_.

- Fetch the latest development versions of the source straight from
  the t3m repository::

    hg clone https://bitbucket.org/t3m/plink
    hg clone https://bitbucket.org/t3m/Spherogram
    hg clone https://bitbucket.org/t3m/SnapPy

- Build the components, from easiest to hardest, and then test::

    cd plink
    python setup.py pip_install
    cd ../Spherogram
    python setup.py pip_install
    cd ../SnapPy
    python setup.py pip_install
    cd ..
    python -m SnapPy.test

- To run the app, you can just do::

    python -m snappy.app

- To build the binary installer, you need PyInstaller::
  
    python -m pip install pyinstaller

  To build the binary installer do::

    cd windows_exe
    python make.py

  You will need to close the SnapPy window that pops up here to
  complete the build process. 

- Useful tips for those coming from Unix.  In MSYS2, your home
  directory is really something like::

    c:\msys2\home\Nathan Dunfield

  whereas your Windows 10 home directory is::

    c:\Users\Nathan Dunfield

  It is handy to have symbolic links from your MSYS2 home directory to
  the Downloads and Desktop folders on the Windows side.  See::
  
    http://www.howtogeek.com/howto/16226/

  for a discussion, but basically you start a "Command Prompt" as
  Adminstrator and do::

    cd "C:\msys64\home\Nathan Dunfield"
    mklink /D Desktop "C:\Users\Nathan Dunfield\Desktop"
    mklink /D Downloads "C:\Users\Nathan Dunfield\Downloads"


- To build CyPari, first install the 32-bit gcc compiler::

    pacman -S mingw-w64-i686-gcc

   Then open a *MinGW32 terminal window*, which is **different** than a
   MSYS2 terminal, and can be started via `c:\msys64\mingw32.exe`.
   This will put the 32-bit gcc in your path and set the correct
   "uname".  Now do::

     hg clone https://bitbucket.org/t3m/CyPari
     cd CyPari
     sh build_pari.sh
     python setup.py build --compiler=mingw32
     python setup.py install
     python -m cypari.test   # About 30 tests will fail.

   Warning: CyPari will not build if there are spaces in the path to
   the CyPari directory.  
