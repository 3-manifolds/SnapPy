Development Basics
==================

Submitting patches
------------------

The source code for SnapPy and its various parts are `hosted on GitHub
<https://github.com/3-manifolds>`_ as `Git repositories
<https://git-scm.com/>`_.  To contribute a patch, create a free
GitHub account, fork the appropriate repository, and then send us a
pull request, as described in `here
<https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request>`_.


macOS
-----

Here is how to get a clean development setup under macOS.

- Install the latest Python 3.13 from Python.org using the `Mac
  Installer Disk Image <http://www.python.org/download/>`_.  Set your
  path so that "python3" is::
      
    /Library/Frameworks/Python.framework/Versions/3.13/bin/python3

- Use `pip <https://pip.pypa.io/en/latest/index.html>`_ to install the
  following packages::

    python3 -m pip install --upgrade setuptools virtualenv wheel pip
    python3 -m pip install cython        # Used for Python-C interfacing
    python3 -m pip install sphinx        # For building the documentation
    python3 -m pip install ipython       # Improved Python shell
    python3 -m pip install py2app        # For making app bundles
    python3 -m pip install pyx FXrays low_index

- Get the source code from the repository, using Git. For
  example you can install Git via its `package installer
  <https://www.git-scm.org/>`_.  Now do::

    git clone https://github.com/3-manifolds/plink.git
    git clone https://github.com/3-manifolds/spherogram.git
    git clone https://github.com/3-manifolds/snappy.git

- Test the stand-alone link editor::

    cd plink
    python3 setup.py pip_install
    python3 -m plink.app   # Link editor appears!

  This last command runs the script ``plink/app.py``; the real code for
  the link editor is in ``plink/__init__.py``.

- Build and install Spherogram::

    cd ../spherogram
    python3 setup.py pip_install
    python3 setup.py test

- Now build SnapPy itself.  One builds it twice to generate the
  documentation, much of which is extracted from the installed module::

    cd ../snappy
    python3 setup.py pip_install
    python3 setup.py test   # Run the tests; pretty picture should appear.
    python3 -m snappy.app   #SnapPy starts!

  To build the clickable app, just do the following::

    cd mac_osx_app
    python3 release.py --manual

  though for general development purposes ``python3 -m snappy.app`` is
  usually the way to go.
    
Some major parts of the SnapPy codebase are:

- ``SnapPy.pyx``, ``SnapPycore.pxi``, ``SnapPy.pxi``: The Cython interface
  to the SnapPea kernel
- ``opengl/CyOpenGL*.pyx``: The Cython interface to OpenGL*
- ``snappy/app.py``: The core GUI code
- ``snappy/polyviewer.py``: The GUI code for Dirichlet domains
- ``snappy/horoviewer.py``: The GUI code for horoball pictures
- ``snappy/database.py``: Interacts with the sqlite3 manifold database


Windows
-------

These instructions have been tested on Windows 10, and describe
setting up a development environment using the (free) MSVC
compiler. To build the CyPari subcomponent, which few will want or
need to do, one must install additional tools.

- Install `Python 3.13.1
  <https://www.python.org/downloads/windows/>`_, specifically the
  default 64-bit version (the file name will end in ``amd64.exe``).
  These instructions assume it has been installed in the directory
  ``C:\Python313`` which is not the default.

- With Python 3.13.1, you need the MSVC command line tools.  You can
  get them by using Microsoft's free `Build Tools for Visual Studio 2019
  <https://visualstudio.microsoft.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=16>`_
  installer and selecting the "C++ command line tools" and installing
  the following subpackages: "MSVC v142 - VS 2019 C++ build tools",
  "Testing tools core features", "C++/CLI support for v142", and
  "Windows 10 SDK (most recent version)".
  
  If you would like to make your own installer you will also need
  `Inno Setup <http://www.jrsoftware.org/isdl.php>`_, specifically the
  unicode version; tested with version 5.5.9.

- Install `MSYS2 <http://msys2.github.io>`_ as appropriate for your
  version Windows.  Follow the instructions on the webpage to update
  everything to the very latest MSYS2 (``pacman -Sy pacman; pacman
  -Syu; pacman -Su`` etc.).  You should also install the ucrt toolchain::

    pacman -S mingw-w64-ucrt-x86_64-toolchain

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

    PATH=/c/Python313:/c/Python313/Scripts:$PATH
    PATH=$PATH:'/c/Program Files (x86)/Inno Setup 5'
    export PATH
    alias emacs="/mingw64/bin/emacs"
    winpty bash; exit

  For example, do::

    nano ~/.bash_profile

- Make sure you have the right version of Python in your path by
  typing::

    python --version

  You should see something like ``Python 3.13.1``.

- Use pip to install some basic tools::
  
    python -m pip install --upgrade pip setuptools wheel  # Upgrades pip to the current version.
    python -m pip install pyreadline sphinx cython cypari pyx FXrays low_index

- Fetch the latest development versions of the source straight from
  the t3m repository::

    git clone https://github.com/3-manifolds/plink.git
    git clone https://github.com/3-manifolds/spherogram.git
    git clone https://github.com/3-manifolds/snappy.git

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
  the Downloads and Desktop folders on the Windows side.  `See this
  discussion <http://www.howtogeek.com/howto/16226/>`_, but basically
  you start a "Command Prompt" as Administrator and do::

    cd "C:\msys64\home\Nathan Dunfield"
    mklink /D Desktop "C:\Users\Nathan Dunfield\Desktop"
    mklink /D Downloads "C:\Users\Nathan Dunfield\Downloads"


Linux
-----

Things you'll need:

- Python 3 with Tkinter: You'll need to have `Python
  <http://python.org>`_ (version 3.9 or newer) and `Tk
  <http://tcl.tk>`_ (at least version 8.5) with `Tkinter
  <http://wiki.python.org/moin/TkInter>`_ to connect them, including
  the header files.  For instance, on Debian or Ubuntu, install the
  packages "python3-tk", "python3-pip", and "python3-dev". On Fedora,
  you'll want e.g. "python3-tkinter", "python3-pip", and
  "python3-devel", and "python3-wheel".

- Test that Python is in order by installing PLink from source::

      python3 -m pip install --user plink
      python3 -m plink.app  # Should start the link editor!

.. _openglmesa:

- Support for OpenGL (3D graphics): This is built in on OS X and the
  most installations of Fedora and Ubuntu.  But you'll need the `MESA
  <http://www.mesa3d.org/>`_ header files "gl.h" and "glu.h" to compile
  SnapPy.  On Debian and Ubuntu, install "libglu1-mesa-dev"; On Fedora install
  "mesa-libGLU-devel".

- `Cython <http://cython.org>`_, which you can install via::

    python3 -m pip install --user cython

- The gcc C++ compiler, g++.

- Fetch the latest development versions of the source straight from
  the repository::

    git clone https://github.com/3-manifolds/PLink.git
    git clone https://github.com/3-manifolds/Spherogram.git
    git clone https://github.com/3-manifolds/Snappy.git

- Build the components, from easiest to hardest, and then test::

    cd PLink
    python setup.py pip_install
    cd ../Spherogram
    python setup.py pip_install
    cd ../SnapPy
    python setup.py pip_install
    cd ..
    python -m SnapPy.test
