Windows
-------------------------------------------------

These instructions have been tested on Windows 7 and 10, and describe
setting up a development environment using the (free) MSVC
compiler. To build the CyPari subcomponent, which few will want or
need to do, one must install additional tools as described at the end.

- Install `Python 2.7 <https://www.python.org/downloads/windows/>`_,
  specifically the 32 bit version (Windows x86 not Windows x86-64).
  Tested with version 2.7.11.

- Install `Python-specific free version of Microsoft Visual C++
  <http://www.microsoft.com/en-us/download/details.aspx?id=44266>`_.
  If you would like to make your own installer you will also need
  `Inno Setup <http://www.jrsoftware.org/isdl.php>`_; tested with
  version 5.5.9.

- Install whichever version of `MSYS2 <http://msys2.github.io>`_ is
  appropriate for your install of Windows.  For concreteness the rest
  of these instructions assume a 64bit version of windows.

- Make a shortcut to `c:\msys64\msys2.exe` as you will be using it all
  the time.

- Install some additional packages::

    pacman -S ssh winpty tar make

- Install your favorite text editor, for example via::

    pacman -S  mingw-w64-x86_64-emacs

- Make it so that MinGW, Python, and Inno Setup are all in
  your PATH, as well as work around some stupid bug, 
  by adding the below lines to the file "~/.bash_profile"::

    alias python="winpty python"
    PATH=/c/Python27:/c/Python27/Scripts:/c/msys64/usr/bin:/c/msys64/mingw64/bin:$PATH
    PATH=$PATH:'/c/Program\ Files\ \(x86\)/Inno\ Setup\ 5
    export PATH
    
  For example, assuming you installed emacs before do::

    /c/msys64/mingw64/bin/emacs ~/.bash_profile

- Python 2.7.9 and newer include `pip
  <https://pip.pypa.io/en/latest/index.html>`_ so let's use it
  to install the needed packages.::
  
    pip install --upgrade pip setuptools     # Upgrades pip to the current version.
    pip install pyreadline sphinx ipython cython
    pip install pyx==0.12.1
    pip install "pyinstaller<3.2"   # There's a bug in 3.2, later versions should be OK.
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

  but to build the binary installer do::

    cd windows_exe
    python make.py

  You will need to close the SnapPy window that pops up here to
  complete the build process. 

- To build CyPari, first install the 32bit gcc compiler::

    pacman -S mingw-w64-i686-gcc

   Then open a MinGW32 terminal window, which is *different* than a
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
