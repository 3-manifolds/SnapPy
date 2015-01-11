.. Installing SnapPy

Installing SnapPy
======================================================

Here are detailed instructions on how to get SnapPy working on a
variety of platforms.  The current version is |release| which was released
on |release_date|.  

Mac OS X
---------------

Simply download `SnapPy.dmg
<http://snappy.computop.org/get/SnapPy.dmg>`_ and copy SnapPy.app
to the Applications folder.  Double-click to start it, just like any
other application.  Works with 10.5 and newer on Macs with Intel
processors.  Users of 10.4 or PPC processors should install `SnapPy-1.4.dmg
<http://snappy.computop.org/get/SnapPy-1.4.dmg>`_ instead.

Windows
-------------------

Simply download and run
`InstallSnapPy.exe. <http://snappy.computop.org/get/InstallSnapPy.exe>`_

NOTE: The Windows version of SnapPy depends on the Microsoft Distributable
Visual C++ Runtime.  If you receive an error message saying
"This application has failed to start because MSVCR90.DLL was not found" or "This application failed to start because the application configuration is incorrect" try downloading and installing `vcredist_x86.exe
<http://www.microsoft.com/downloads/details.aspx?FamilyID=9b2da534-3e03-4391-8a4d-074b9f2bc1bf&displaylang=en>`_ from Microsoft.

If you are running Windows 7 and the program works except for the 3D
graphics features, then you are likely missing "msvcr71.dll".

Linux
--------------------

Here are short recipes which work on many Linux systems, with both
32-bit and 64-bit kernels supported. These instructions assume you
have system administrator (superuser) privileges; if not, you can
install SnapPy into a `virtual environment`_ *assuming* the needed
packages are installed.  For other systems, try the one closest to
yours below, and if that fails, follow the instructions for `generic
Unix`_ in the next section.

+ **Fedora:** Tested on versions 8-10, 14 (Werewolf-Sulfer-Cambridge, Laughlin)::

    sudo yum install tkinter python-setuptools-devel 
    sudo python -m easy_install -U -f http://snappy.computop.org/get snappy

  Note: For this to work, you need to set the SELinux Enforcement mode
  to Permissive or lower.

+ **Ubuntu/Debian:** Tested on Ubuntu 8.04 (Hardy Heron), 8.10 (Intrepid Ibex), 9.04 (Jaunty Jackalope), 9.10 (Karmic Koala), 10.10 (Maverick Meerkat)::

    sudo apt-get install python-tk python-setuptools    
    sudo python -m easy_install -U -f http://snappy.computop.org/get snappy

+ **PCLinuxOS:** Not actually tested, but should work::

    sudo apt-get install tkinter python-setuptools
    sudo python -m easy_install -U -f http://snappy.computop.org/get snappy

Once you have it installed, do::

  python -m snappy.app

You may get a message about creating a ".ipython" directory; this is
normal, just hit return to continue.  There should also now be a
command "SnapPy" which does the same thing.

Note: You need to have Python 2.6 or 2.7 to install SnapPy 1.4.0 or
newer; if instead you have Python 2.5 the above instructions will
install 1.3.12 instead.


Generic Unix
----------------------------------------------------------

If you use a Unix other that OS X or Linux, or if the prebuilt
packages don't work for you, you'll need to build SnapPy from source.
Here are some detailed instructions.

Things you'll need:

- Python with Tkinter: You'll need to have `Python <http://python.org>`_
  (version 2.6-2.7) and `Tk <http://tcl.tk>`_ (at least version 8.4)
  with `Tkinter <http://wiki.python.org/moin/TkInter>`_ to
  connect them, including the header files.  For instance, on Debian
  or Ubuntu, install the packages "python-tk" and "python-dev". On
  Fedora, you'll want "tkinter" and "python-devel". In addition, you'll
  need

  - `Setuptools <http://pypi.python.org/pypi/distribute>`_, which is
    typically packaged as "python-setuptools" (Ubuntu/Debian),
    "python-setuptools-devel" (Fedora), or can be installed via::

      curl -O http://peak.telecommunity.com/dist/ez_setup.py
      sudo python ez_setup.py  

    Test that Python is in order by installing PLink from source::

      python -m easy_install -f http://t3m.computop.org/plink plink
      plink   # Should start the link editor!

.. _openglmesa:

- Support for OpenGL (3D graphics): This is built in on OS X and the
  most installations of Fedora and Ubuntu.  But you'll need the `MESA
  <http://www.mesa3d.org/>`_ header files "gl.h" and "glu.h" to compile
  SnapPy.  On Debian and Ubuntu, install "libglu1-mesa-dev"; On Fedora install
  "mesa-libGLU-devel".

- `Cython <http://cython.org>`_, which you can install via::

    sudo python -m easy_install cython

- `Sphinx <http://sphinx.pocoo.org/>`_, which you can install via::

    sudo python -m easy_install sphinx

- The gcc C++ compiler, g++, which is not installed by default on some
  systems, e.g. Ubuntu 11.10.

- `CyPari <http://www.math.uic.edu/t3m/>`_: a stand-alone version of
  `Sage's <http://sagemath.org>`_ Python interface to the
  `PARI <http://pari.math.u-bordeaux.fr/PARI>`_ number theory library.

Now download the `source code`_ listed below, for instance::

    curl -L -O http://snappy.computop.org/get/SnapPy.tar.gz
    tar xfz SnapPy.tar.gz; cd SnapPy

There is one more dependency that need to be dealt with:

- `Togl <http://togl.sf.net>`_: a 3d widget for Tk. For OS X and
  Linux, there are pre-built binaries of this in the snappy
  subdirectory, e.g. snappy/linux2-tk8.4.  For Linux these are built for
  both 32-bit and 64-bit kernels, and should work on most systems.  If
  they don't, you'll need to edit or follow "build_togl.sh" to build
  Togl directly into "snappy/linux2-tk*" (32-bit kernel) or
  "snappy/linux2-x86_64-tk*" (64-bit kernel), where "*" is the version
  of Tk you are using.
  
Finally, compile and install the SnapPy module (which will install
certain other dependencies) and test::

  sudo python setup.py install
  sudo python setup.py build_docs install
  cd /tmp; python -m snappy.app

You may get a message about creating a ".ipython" directory; this is
normal, just hit return to continue.  There should also now be a
command "SnapPy" which does the same thing as "python -m snappy.app".

Python Modules for Macintosh or Windows
---------------------------------------

If you write Python programs on a Macintosh or Windows system, you may
wish to install SnapPy as a Python module into your own copy of Python
2.6 (Mac only) or 2.7 (both platforms). After installing Python and
`setuptools <http://pypi.python.org/pypi/distribute>`_, you may
install a SnapPy module from your Terminal application or Command
Prompt with the command::

    python -m easy_install -U -f http://snappy.computop.org/get snappy

OS X notes: For best results, use a Python downloaded from `Python.org
<http://python.org>`_ and not the one provided by Apple.  You need at
least 10.5 and an Intel processor to use the latest versions of these
precompiled modules; if you have an older system, you will get version
1.4.* instead.


Virtual Environment
-----------------------------------

All of the above instructions assume that you want to install SnapPy
globally, in the main Python site-packages directory.  You can also
create a Python `virtual environment <http://www.virtualenv.org/>`_
and install SnapPy into it.  For example, to install SnapPy into
"mypy/bin" do::

   #Download needed files, could also use any webbrowser here.
   curl -O https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.2.tar.gz
   tar xfz virtualenv-*.tar.gz
   # Create a virtual environment in new directory "mypy" 
   python virtualenv-*/virtualenv.py mypy 
   # Install and run SnapPy!
   mypy/bin/easy_install -U -f http://snappy.computop.org/get  snappy
   mypy/bin/SnapPy

Sage
----

SnapPy has some special features when used within `Sage
<http://sagemath.org>`_, the universal mathematics software based on
Python.   You can install it as a Sage optional package via

.. parsed-literal::

   sage -i \http://snappy.computop.org/get/snappy-|release|.spkg

and as of December 2013 the version of Sage on the `Sagemath Cloud
<https://cloud.sagemath.com/>`_ has SnapPy preinstalled.  

If it has trouble when compiling CyOpenGL, you are probably missing
the `"gl.h" and "glu.h" headers <installing.html#openglmesa>`_.  The graphical
features may or may not work, depending on how Tkinter was configured
within Sage.  If you are using Sage 5.11 or newer, the graphics
features may seem to "hang" when you try to start them.  If this
happens, type "%gui none" at the Sage prompt; please note that doing so
will break Sage's "attach" feature.

Source code
-----------------------------------

The complete source code for all platforms: `SnapPy.tar.gz <http://snappy.computop.org/get/SnapPy.tar.gz>`_   

You can also get it straight from the `Mercurial
<http://www.selenic.com/mercurial>`_ repository::

  hg clone static-http://math.uic.edu/t3m/hg/SnapPy




