Windows
-------

These instructions have been tested on Windows 10, and describe
setting up a development environment using the (free) MSVC
compiler. To build the CyPari subcomponent, which few will want or
need to do, one must install additional tools as described at the end.

- Install `Python 2.7 <https://www.python.org/downloads/windows/>`_,
  specifically the 32 bit version (Windows x86 not Windows x86-64).
  Tested with version 2.7.13.

- Install `Python-specific free version of Microsoft Visual C++
  <http://www.microsoft.com/en-us/download/details.aspx?id=44266>`_.
  If you would like to make your own installer you will also need
  `Inno Setup <http://www.jrsoftware.org/isdl.php>`_; tested with
  version 5.5.9.

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

    pacman -S git make nano openssh perl tar unzip wget winpty

- Install your favorite text editor, for example you can install Emacs
  via::

    pacman -S  mingw-w64-x86_64-emacs

- Make it so that MinGW, Python, and Inno Setup are all in your PATH,
  as well as work around some stupid bug, by making the end of your
  "~/.bash_profile" file to read::

    PATH=/c/Python27:/c/Python27/Scripts:$PATH
    PATH=$PATH:'/c/Program Files (x86)/Inno Setup 5'
    export PATH
    alias python="winpty python"
    alias ipython="winpty ipython"
    alias pip="winpty pip"
    alias emacs="/mingw64/bin/emacs"
    alias dumpbin="/c/Users/Nathan\ Dunfield/AppData/Local/Programs/Common/Microsoft/Visual C++ for Python/9.0/VC/bin/dumpbin.exe"

  For example, do::

    nano ~/.bash_profile

- Python 2.7.9 and newer include `pip
  <https://pip.pypa.io/en/latest/index.html>`_ so let's use it
  to install the needed packages.::
  
    pip install --upgrade pip setuptools     # Upgrades pip to the current version.
    pip install pyreadline sphinx cython cypari
    pip install pyx==0.12.1
    pip install mercurial   # Installs "hg", used in next step

- Fetch the latest development versions of the source straight from
  the t3m repository::

    hg clone https://bitbucket.org/t3m/plink
    hg clone https://bitbucket.org/t3m/Spherogram
    hg clone https://bitbucket.org/t3m/SnapPy

- Build the components, from easiest to hardest, and then test::

    cd plink
    python setup.py install
    cd ../Spherogram
    python setup.py install
    cd ../SnapPy
    python setup.py install
    cd ..
    python -m SnapPy.test

- To run the app, you can just do::

    python -m snappy.app

- To build the binary installer, you need PyInstaller, but you need a
  special copy::

    pip install "pyinstaller<3.2"

  but to build the binary installer do::

    cd windows_exe
    python make.py

  You will need to close the SnapPy window that pops up here to
  complete the build process. 

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
