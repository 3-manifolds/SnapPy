.. Installing SnapPy

Installing SnapPy
=================

Here are detailed instructions on how to get SnapPy working on a
variety of platforms.  The current version is |release| which was released
on |release_date|.  If you encounter problems installing SnapPy,
`please let us know <bugs.html>`_. 

macOS/OS X
----------

Simply download `SnapPy.dmg
<https://bitbucket.org/t3m/snappy/downloads/SnapPy.dmg>`_ and copy
SnapPy.app to the Applications folder.  Double-click to start it, just
like any other application.  Works with OS X versions 10.6 and newer.
Users of earlier versions of OS X should install `SnapPy-1.4.dmg
<http://t3m.computop.org/SnapPy-nest/SnapPy-1.4.dmg>`_ instead.  We
also offer a version based on Python 3 rather than the default
Python 2: `SnapPy-Python3.dmg
<https://bitbucket.org/t3m/snappy/downloads/SnapPy-Python3.dmg>`_.

Windows
-------

Simply download and run
`InstallSnapPy.exe <https://bitbucket.org/t3m/snappy/downloads/InstallSnapPy.exe>`_.
We also offer a version based on Python 3 rather than the default
Python 2: `InstallSnapPy-Python3.exe
<https://bitbucket.org/t3m/snappy/downloads/InstallSnapPy-Python3.exe>`_.


Linux
-----

Here are short recipes which work on most Linux systems, specifically
those that run a 64-bit kernel and has Python 3.4 or newer. These
instructions assume you have system administrator (superuser)
privileges to install software packages from your Linux distribution
but want to install SnapPy (and its dependencies) in just your user
directory, specifically ``~/.local``.  For other Linux systems, try
the one closest to yours below, and if that fails, follow the
instructions for `generic Unix`_ in the next section.

+ **Ubuntu/Debian/Mint**:: Tested on Ubuntu 16.04 and 18.04 and Debian 

    sudo apt-get install python3-tk python3-pip
    # Note no "sudo" on the next one!
    python3 -m pip install --upgrade --user snappy

  To run the app, do::

    ~/.local/bin/SnapPy
    
+ **Red Hat Enterprise Linux/CentOS/SciLinux**:: These instructions
  are for version 7 or later, and you need to have the `EPEL packages
  available
  <https://fedoraproject.org/wiki/EPEL#How_can_I_use_these_extra_packages.3F>`_.
  For CentOS and SciLinux, you can access EPEL packages by doing::

    sudo yum install epel-release

  Now install via::
    
    sudo yum install python36-tkinter python36-pip
    # Note no "sudo" on the next one!
    python36 -m pip install --upgrade --user snappy

  To run the app, do::

    ~/.local/bin/SnapPy

+ **Arch/Manjaro**::

+ **openSUSE**::

+ **PCLinuxOS:** Untested, but try the instructions for Ubuntu.  

Once you have installed SnapPy, do::

  python -m snappy.app

You may get a message about creating a ".ipython" directory; this is
normal, just hit return to continue.  There should also now be a
command "SnapPy" which does the same thing.



Generic Unix
------------

If you use a Unix other than OS X or Linux, or if the prebuilt
packages don't work for you, you'll need to build SnapPy from source.
Here are some detailed instructions.

Things you'll need:

- Python with Tkinter: You'll need to have `Python <http://python.org>`_
  (version 2.6-2.7) and `Tk <http://tcl.tk>`_ (at least version 8.4)
  with `Tkinter <http://wiki.python.org/moin/TkInter>`_ to
  connect them, including the header files.  For instance, on Debian
  or Ubuntu, install the packages "python-tk" and "python-dev". On
  Fedora, you'll want "tkinter" and "python-devel". In addition, you'll
  need `setuptools <https://pypi.python.org/pypi/setuptools>`_, which is
  typically packaged as "python-setuptools".

- Test that Python is in order by installing PLink from source::

      python -m easy_install plink
      plink   # Should start the link editor!

.. _openglmesa:

- Support for OpenGL (3D graphics): This is built in on OS X and the
  most installations of Fedora and Ubuntu.  But you'll need the `MESA
  <http://www.mesa3d.org/>`_ header files "gl.h" and "glu.h" to compile
  SnapPy.  On Debian and Ubuntu, install "libglu1-mesa-dev"; On Fedora install
  "mesa-libGLU-devel".

- `Cython <http://cython.org>`_, which you can install via::

    sudo python -m easy_install cython

- The gcc C++ compiler, g++, which is not installed by default on some
  systems, e.g. Ubuntu 11.10.

- `CyPari <https://pypi.python.org/pypi/cypari/>`_: a stand-alone version of
  `Sage's <http://sagemath.org>`_ Python interface to the
  `PARI <http://pari.math.u-bordeaux.fr/PARI>`_ number theory library.

Now download the `source code`_ listed below, for instance

.. parsed-literal::
   
   curl -L -O |tarball|  
   tar xfz |tarball|; cd SnapPy

There is one more dependency that may need to be dealt with:

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
  python -m snappy.app

You may get a message about creating a ".ipython" directory; this is
normal, just hit return to continue.  There should also now be a
command "SnapPy" which does the same thing as "python -m snappy.app".

Python Modules for Macintosh or Windows
---------------------------------------

If you write Python programs on a Macintosh or Windows system, you may
wish to install SnapPy as a Python module into your own copy of Python
2.7.  (For best results on macOS, use a Python downloaded from
`Python.org <http://python.org>`_ and not the one provided by Apple.)
After installing Python, you may install a SnapPy module from your
Terminal application or Command Prompt with the commands::

    python -m pip install --upgrade pip setuptools
    python -m pip install --upgrade --upgrade-strategy only-if-needed snappy

If your Python lacks the pip module, `get it here
<https://pip.pypa.io/en/stable/installing/>`_.


Virtual Environment
-------------------

All of the above instructions assume that you want to install SnapPy
globally, in the main Python site-packages directory.  You can also
create a `Python virtual environment <http://www.virtualenv.org/>`_
and install SnapPy into it.  For example, to install SnapPy into
"mypy/bin" do::

   # Create a virtual environment in new directory "mypy" 
   python -m virtualenv mypy 
   # Install and run SnapPy!
   mypy/bin/easy_install snappy
   mypy/bin/SnapPy

SageMath
--------

SnapPy has some special features when used within `SageMath
<http://sagemath.org>`_, the universal mathematics software based on
Python.   You can install it as a Sage optional package via the
following if using Sage 6.4 or newer::

  sage -pip install snappy

If you are on macOS and it complains about not having SSL, TLS, or
something related to a certificate missing, you likely have the
problem `described here
<https://groups.google.com/d/msg/sage-devel/h974Gv6kOtg/XDJj9ByiBgAJ>`_
so try `this approach
<https://groups.google.com/d/msg/sage-devel/h974Gv6kOtg/Fq49Qo3vBgAJ>`_
If you encounter other problems, on any platform, try::

  sage -pip install --no-binary :all: snappy

For Sage 6.3 or older do::
  
  sage -python -m easy_install snappy

Alternatively, SageMath on `CoCalc <https://cocalc.com/>`_ (formerly
the SageMathCloud) also has SnapPy preinstalled, and the graphics
features even work via the `X11 interface
<http://blog.sagemath.com/cocalc/2018/11/05/x11.html>`_, see the
bottom of that page for more.

If you previously installed SnapPy into SageMath and want to upgrade
SnapPy to the latest version, do::

  sage -pip install --upgrade --no-deps snappy_manifolds plink spherogram FXrays decorator snappy

or::

  sage -python -m easy_install -U snappy

as appropriate.

If it has trouble when compiling CyOpenGL, you are probably missing
the `"gl.h" headers <installing.html#openglmesa>`_.  The graphical
features may or may not work, depending on how Tkinter was configured
within Sage, and may seem to "hang" when you try to start them.  To
deal with the latter issue on Sage 5.11 or later, type "%gui tk" at
the Sage prompt; please note that doing so may break Sage's "attach"
feature.

Source code
-----------

The complete source code for all platforms: |tarball|_

You can also browse our `source code repository
<https://bitbucket.org/t3m/snappy>`_ or clone it using `Mercurial <http://mercurial-scm.org/>`_ via::

  hg clone https://bitbucket.org/t3m/snappy

Python 3
--------

We now fully support using SnapPy with Python 3!  Currently, binaries
are provided for Python 3.4, 3.5, and 3.6 on macOS, Linux, and
Windows.  We offer stand-alone application for macOS (`SnapPy-Python3.dmg
<https://bitbucket.org/t3m/snappy/downloads/SnapPy-Python3.dmg>`_) and
Windows (`InstallSnapPy-Python3.exe
<https://bitbucket.org/t3m/snappy/downloads/InstallSnapPy-Python3.exe>`_).
You can also install the Python modules into your Python via the
following, with Linux users needing to add ``sudo`` at the start of
each line::

  python3 -m pip install --upgrade pip setuptools
  python3 -m pip install --upgrade --upgrade-strategy only-if-needed snappy
  python3 -m snappy.app
