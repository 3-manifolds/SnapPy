.. Installing SnapPy

Installing SnapPy
=================

Here are detailed instructions on how to get SnapPy working on a
variety of platforms.  The current version is |release| which was released
on |release_date|.  If you encounter problems installing SnapPy,
`please let us know <bugs.html>`_. 

macOS
-----

Simply download `SnapPy.dmg
<https://bitbucket.org/t3m/snappy/downloads/SnapPy.dmg>`_ and copy
SnapPy.app to the Applications folder.  Double-click to start it, just
like any other application.  Works with macOS/OS X versions 10.6 and
newer.  Users of earlier versions of OS X should install
`SnapPy-1.4.dmg <https://t3m.computop.org/SnapPy-nest/SnapPy-1.4.dmg>`_
instead.  We also offer a version based on Python 2 rather than the
default Python 3: `SnapPy-Python2.dmg
<https://bitbucket.org/t3m/snappy/downloads/SnapPy-Python2.dmg>`_;
this was the default download in SnapPy 2.6.1 and earlier.

Windows
-------

Simply download and run
`InstallSnapPy.exe <https://bitbucket.org/t3m/snappy/downloads/InstallSnapPy.exe>`_.
We also offer a version based on Python 2 rather than the default
Python 3: `InstallSnapPy-Python2.exe
<https://bitbucket.org/t3m/snappy/downloads/InstallSnapPy-Python2.exe>`_.


Linux
-----

Here are short recipes which work on most Linux systems, specifically
those that run a 64-bit kernel and have Python 3.4 or newer. These
instructions assume you have system administrator (superuser)
privileges to install software packages from your Linux distribution
but want to install SnapPy (and its various Python dependencies) just
in your own user directory, specifically ``~/.local``.  For other
Linux systems, try the one closest to yours below, and if that fails,
follow the instructions for `generic Unix`_.

+ **Ubuntu/Debian/Mint**: Tested on Ubuntu 16.04 and 18.04 and Debian::

    sudo apt-get install python3-tk python3-pip
    # Note no "sudo" on the next one!
    python3 -m pip install --upgrade --user snappy

+ **Fedora**: Tested on Fedora 30::

    sudo yum install python3-tkinter python3-pip
    # Note no "sudo" on the next one!
    python3 -m pip install --upgrade --user snappy
    
+ **Red Hat Enterprise Linux/CentOS/SciLinux**: These instructions
  are for version 7 or later, and you need to have the `EPEL packages
  available
  <https://fedoraproject.org/wiki/EPEL#How_can_I_use_these_extra_packages.3F>`_.
  For CentOS and SciLinux, you can access EPEL packages by doing::

    sudo yum install epel-release

  Now install via::
    
    sudo yum install python36-tkinter python36-pip
    # Note no "sudo" on the next one!
    python36 -m pip install --upgrade --user snappy

+ **Arch/Manjaro**: Install via::

    sudo pacman -Sy python-pip tk
    # Note no "sudo" on the next one!
    python -m pip install --upgrade --user snappy

+ **openSUSE**: Install via::

    sudo zypper install -y python3-tk python3-pip
    # Note no "sudo" on the next one!
    python3 -m pip install --upgrade --user snappy

Once you have installed SnapPy, do the following to start it::

    ~/.local/bin/SnapPy

You may get a message about creating a ".ipython" directory; this is
normal, just hit return to continue.  There should also now be a
command "SnapPy" which does the same thing.  To make it so that you
can start SnapPy with just the command ``SnapPy``, make sure
``~/.local/bin`` is in `in your path
<https://opensource.com/article/17/6/set-path-linux>`_.


Python Modules for Macintosh or Windows
---------------------------------------

If you write Python programs on a Macintosh or Windows system, you may
wish to install SnapPy as a Python module into your own copy of
Python.  We support Python 3.4 and up as well as the legacy Python
2.7.  (On macOS, use a Python downloaded from `Python.org
<http://python.org>`_ and not the one provided by Apple.)  After
installing Python, you may install a SnapPy module from your Terminal
application or Command Prompt with the commands::

    python3 -m pip install --upgrade --user snappy

If you use Python 2 rather than Python 3, replace ``python3`` with
``python`` in the above.  If your Python lacks the pip module, `get it
here <https://pip.pypa.io/en/stable/installing/>`_.


SageMath
--------

SnapPy has some special features when used within `SageMath
<http://sagemath.org>`_, the universal mathematics software based on
Python.  This section describes how to install SnapPy into your
existing copy of SageMath, but you may find it easier to use the
`kitchen sink`_ approach instead.  You can install it as a Sage
optional package via the following if using Sage 6.4 or newer::

  sage -pip install snappy

If you are on macOS and it complains about not having SSL, TLS, or
something related to a certificate missing, you likely have the
problem `described here
<https://groups.google.com/d/msg/sage-devel/h974Gv6kOtg/XDJj9ByiBgAJ>`_
so try `this approach
<https://groups.google.com/d/msg/sage-devel/h974Gv6kOtg/Fq49Qo3vBgAJ>`_
If you encounter other problems, on any platform, try::

  sage -pip install --no-binary :all: snappy

Alternatively, SageMath on `CoCalc <https://cocalc.com/>`_ (formerly
the SageMathCloud) also has SnapPy preinstalled, and the graphics
features even work via the `X11 interface
<http://blog.sagemath.com/cocalc/2018/11/05/x11.html>`_, see the
bottom of that page for more.

If you previously installed SnapPy into SageMath and want to upgrade
SnapPy to the latest version, do::

  sage -pip install --upgrade --no-deps snappy_manifolds plink spherogram FXrays decorator snappy

If it has trouble when compiling CyOpenGL, you are probably missing
the `"gl.h" headers <installing.html#openglmesa>`_.  The graphical
features may or may not work, depending on how Tkinter was configured
within Sage, and may seem to "hang" when you try to start them.  To
deal with the latter issue type "%gui tk" at the Sage prompt; please
note that doing so may break Sage's "attach" feature.


Kitchen sink
------------

SnapPy gains extra features when used in `SageMath`_ and one can use
Sage's Python to interact not just with SnapPy but a range of other
computational tools in low-dimensional topology including
`Regina <http://regina-normal.github.io/>`_,
`snap <http://snap-pari.sourceforge.net>`_,
`heegaard <https://bitbucket.org/t3m/heegaard>`_,
`gridlink <https://bitbucket.org/t3m/gridlink>`_,
and `flipper <http://flipper.readthedocs.io>`_.
We offer a `prepackaged Docker image
<https://bitbucket.org/t3m/sagedocker>`_ with all of the above tools
and many more; using this is frequently the easiest way to get a
working setup for such multifaceted computations.  For more, watch
`this demonstration <https://icerm.brown.edu/video_archive/?play=1992>`_.

We also offer `conda environments
<https://github.com/unhyperbolic/condaForSnapPy>`_ with SnapPy and
optionally Sage (only on Mac OS and Linux). While it has none of the
other aforementioned tools, it has the advantage that the GUI elements
such as the link editor and the browser can be used directly.

Generic Unix
------------

If you use a Unix other than OS X or Linux, or if the prebuilt
packages don't work for you, you'll need to build SnapPy from source.
Here are some detailed instructions.

Things you'll need:

- Python 3 with Tkinter: You'll need to have `Python
  <http://python.org>`_ (version 3.4 or newer) and `Tk
  <http://tcl.tk>`_ (at least version 8.4) with `Tkinter
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

    python3 -m pip --user cython

- The gcc C++ compiler, g++.

- `CyPari <https://pypi.python.org/pypi/cypari/>`_: a stand-alone version of
  `Sage's <http://sagemath.org>`_ Python interface to the
  `PARI <http://pari.math.u-bordeaux.fr/PARI>`_ number theory
  library.  Usually, you can install this with::

     python3 -m pip install --user cypari

Now download the `source code`_ listed below, for instance

.. parsed-literal::
   
   wget https://pypi.python.org/packages/source/s/snappy/|tarball|  
   tar xfz |tarball|; cd snappy-*

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

  python3 setup.py build
  python3 -m pip install --user .
  python3 -m snappy.test
  python3 -m snappy.app


Source code
-----------

The complete source code for all platforms: |tarball|_

You can also browse our `source code repository
<https://github.com/3-manifolds/SnapPy>`_ or clone it using `git
<https://git-scm.com/>`_ via::

  git clone https://github.com/3-manifolds/SnapPy.git
