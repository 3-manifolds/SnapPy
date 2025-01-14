.. Installing SnapPy

Installing SnapPy
=================

Here are detailed instructions on how to get SnapPy working on a
variety of platforms.  The current version is |release| which was released
on |release_date|.  If you encounter problems installing SnapPy,
:doc:`please let us know <bugs>`.

macOS
-----

Simply download `SnapPy.dmg
<https://github.com/3-manifolds/SnapPy/releases/latest/download/SnapPy.dmg>`_
and copy SnapPy.app to the Applications folder.  Double-click to start
it, just like any other application.  Works with macOS 10.13 and
newer.  Earlier releases `can be found here
<https://github.com/3-manifolds/SnapPy/releases/>`_.

Windows
-------

Simply download and run `InstallSnapPy.exe
<https://github.com/3-manifolds/SnapPy/releases/latest/download/InstallSnapPy.exe>`_.
Earlier releases `can be found here
<https://github.com/3-manifolds/SnapPy/releases/>`_.


Linux app
---------

Starting with SnapPy 3.2, a completely self-contained SnapPy `AppImage
<https://docs.appimage.org/introduction/quickstart.html#ref-quickstart>`_
is available that should work on any Linux system from the last 5
years.  This AppImage contains its own private copy of Python, so if
you plan to use SnapPy in your own Python program skip ahead to
`Python Modules for Linux`_.  Here is the recipe for installing the
AppImage in ``~/bin`` after you have downloaded `SnapPy-x86_64.AppImage
<https://github.com/3-manifolds/SnapPy/releases/latest/download/SnapPy-x86_64.AppImage.>`_::

  mkdir -p ~/bin
  mv ~/Downloads/SnapPy-x86_64.AppImage ~/bin
  chmod +x ~/bin/SnapPy-x86_64.AppImage
  ln -s -f ~/bin/SnapPy-x86_64.AppImage ~/bin/SnapPy
  ~/bin/SnapPy-x86_64.AppImage --install

The last command registers the SnapPy app with your desktop system and
starts SnapPy.  In the future, you can start SnapPy by using the desktop
search tool with Gnome or the main menu with KDE. You can pin the icon
to your Dash or Task Bar for easy access.  From a terminal window,
you can also start the app by typing ``SnapPy`` provided ``~/bin`` is
in your `$PATH <https://opensource.com/article/17/6/set-path-linux>`_.


Python Modules for Linux
------------------------

If you want SnapPy to use the system version of Python, for example to
incorporate SnapPy in your own Python scripts, below are short recipes
for doing this on most common Linux system.  These instructions assume
you have system administrator (superuser) privileges to install
software packages from your Linux distribution.  (If you're not a
superuser, you can still use the `Linux app`_ or try `Conda`_.) For other
Linux systems, try the one closest to yours below, and if that fails,
follow the instructions for `generic Unix`_.

The first step is to install Python and other requirements.

+ **Ubuntu/Debian/Mint/MX Linux/Elementary:** Tested on Ubuntu 24.04::

    sudo apt install python3-pip python3-tk

+ **Fedora**: Tested on Fedora 41::

    sudo dnf install python3-pip python3-tkinter

+ **Arch/Manjaro/EndeavourOS**: Install via::

    sudo pacman -Sy python-pip tk

+ **openSUSE**: For openSUSE Tumbleweed::

    sudo zypper install python3-tk

  For openSUSE Leap, as of verion 15.6 you need ask for a recent
  version of Python or it will give you Python 3.6 which is too old
  for SnapPy::
    
    sudo zypper install python3.12-tk

  You will need to replace ``python3`` by ``python3.12`` in subsequent
  steps.
    
+ **Red Hat Enterprise Linux/CentOS/Rocky Linux/AlmaLinux:**: These instructions
  are for version 8 or later; tested on AlmaLinux 8 and 9::

    sudo dnf install python3.11-pip python3.11-tkinter

  You will need to replace ``python3`` by ``python3.11`` in subsequent
  steps.


Next, you need to install the SnapPy python modules. Ideally, this
would be done in a venv that you use for your SnapPy project.  However
if you are not familiar with virtual environments, or if you have
followed install instructions for previous versions of SnapPy, it may
be easiest and safest to install the new version the same way, namely
by doing a ``user install`` with pip. 

To do a ``user install`` with pip, the first thing to try is::

  # Note no "sudo" below!
  python3 -m pip install --upgrade --user snappy

If you get a long error message that starts::

  error: externally-managed-environment

then you should probably set up a virtual environment and install SnapPy
into it, although an alternative is suggested below.

Here is the `official tutorial
<https://docs.python.org/3/tutorial/venv.html>`_ on using virtual
environments and an `indepth discussion
<https://realpython.com/python-virtual-environments-a-primer/>`_.  A
recipe is::

  python3 -m venv snappy_venv
  # Switch to snappy_venv's Python
  source snappy_venv/bin/activate
  pip install snappy
  # Start SnapPy app!
  SnapPy
  # Return to system Python
  deactivate

If you always want to use the ``snappy_venv`` Python, adjust your `$PATH
<https://opensource.com/article/17/6/set-path-linux>`_ to include ``snappy_venv/bin``.

The alternative way to work around the ``externally-managed-environment``
error is to do the following::

  # Note no "sudo" below!
  python3 -m pip install --upgrade --user --break-system-packages snappy

Despite the scary name, provided you don't use ``sudo``, the flags
``--user --break-system-packages`` will not actually modify the system
packages and will just install ``snappy`` into the subdirectory
``~/.local/share/python3.*/site-packages`` of your home directory,
just as ``--user`` does on more permissive systems.
    
If you want the larger version of HTLinkExteriors that includes the 15
crossing knots (uses 110M of disk space), also install the Python
package ``snappy_15_knots``, for example::

  python3 -m pip install --upgrade --user snappy_15_knots

Once you have installed SnapPy, just run the following command to start
the app::

    ~/.local/bin/SnapPy

So that you can start SnapPy with just the command ``SnapPy``, make
sure ``~/.local/bin`` is in `in your path
<https://opensource.com/article/17/6/set-path-linux>`_.


Python Modules for Macintosh or Windows
---------------------------------------

If you write Python programs on a Macintosh or Windows system, you may
wish to install SnapPy as a Python module into your own copy of
Python.  We support Python 3.6 and up.  (On macOS, use a Python
downloaded from `Python.org <http://python.org>`_ and not the one
provided by Apple.)  After installing Python, you may install a SnapPy
module from your Terminal application or Command Prompt with the
commands::

    python3 -m pip install --upgrade --user snappy

If you want the larger version of HTLinkExteriors that includes the 15
crossing knots (uses 110M of disk space), do::

    python3 -m pip install --upgrade --user snappy_15_knots

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
  sage -pip install snappy_15_knots  # Larger version of HTLinkExteriors

If you are on macOS, we recommend use `this binary
<https://github.com/3-manifolds/Sage_macOS/releases>`_.

Alternatively, SageMath on `CoCalc <https://cocalc.com/>`_ (formerly
the SageMathCloud) also has SnapPy preinstalled, and the graphics
features even work via the `X11 interface
<http://blog.sagemath.com/cocalc/2018/11/05/x11.html>`_, see the
bottom of that page for more.

If you previously installed SnapPy into SageMath and want to upgrade
SnapPy to the latest version, do::

  sage -pip install --upgrade snappy

If it has trouble when compiling CyOpenGL, you are probably missing
the `"gl.h" headers <openglmesa>`.  The graphical features may or may
not work, depending on how Tkinter was configured within Sage, and may
seem to "hang" when you try to start them.  To deal with the latter
issue type "%gui tk" at the Sage prompt; please note that doing so may
break Sage's "attach" feature.


Kitchen sink
------------

SnapPy gains extra features when used in `SageMath`_ and one can use
Sage's Python to interact not just with SnapPy but a range of other
computational tools in low-dimensional topology including
`Regina <http://regina-normal.github.io/>`_,
`snap <http://snap-pari.sourceforge.net>`_,
`heegaard <https://github.com/3-manifolds/heegaard>`_,
`gridlink <https://github.com/3-manifolds/gridlink>`_,
and `flipper <http://flipper.readthedocs.io>`_.
We offer a `prepackaged Docker image
<https://hub.docker.com/r/computop/sage/>`_ with all of the above tools
and many more; using this is frequently the easiest way to get a
working setup for such multifaceted computations.  For more, watch
`this demonstration <https://icerm.brown.edu/video_archive/?play=1992>`_.

We also offer `conda environments
<https://github.com/unhyperbolic/condaForSnapPy>`_ with SnapPy and
optionally Sage (only on Mac OS and Linux). While it has none of the
other aforementioned tools, it has the advantage that the GUI elements
such as the link editor and the browser can be used directly.

Conda
-----

**FILL IN***


Generic Unix
------------

If you use a Unix other than OS X or Linux, or if the prebuilt
packages don't work for you, you'll need to build SnapPy from source.
Here are some detailed instructions.

Things you'll need:

- Python 3 with Tkinter: You'll need to have `Python
  <http://python.org>`_ (version 3.6 or newer) and `Tk
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
  subdirectory, e.g. snappy/linux2-tk8.4.  For Linux these are built
  for 64-bit kernels, and should work on most systems.  If they don't,
  you'll need to edit or follow "build_togl.sh" to build Togl directly.

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
