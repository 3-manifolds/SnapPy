.. Installing SnapPy

Installing SnapPy
=================

Here are detailed instructions on how to get SnapPy working on a
variety of platforms.  The current version is |release| which was released
on |release_date|.  If you encounter any problems installing SnapPy,
:doc:`please let us know <bugs>`.

macOS
-----

Simply download `SnapPy.dmg
<https://github.com/3-manifolds/SnapPy/releases/latest/download/SnapPy.dmg>`_
and copy SnapPy.app to the Applications folder.  Double-click to start
it, just like any other application.  The current version works with macOS 10.14 and
newer and earlier releases `can be found here
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
you plan to use SnapPy in your own Python programs, skip ahead to
`Python Modules for Linux`_.  Here is the recipe for installing the
AppImage in ``~/bin`` after you have downloaded the file
`SnapPy-x86_64.AppImage
<https://github.com/3-manifolds/SnapPy/releases/latest/download/SnapPy-x86_64.AppImage.>`_::

  mkdir -p ~/bin
  mv ~/Downloads/SnapPy-x86_64.AppImage ~/bin
  chmod +x ~/bin/SnapPy-x86_64.AppImage
  ln -s -f ~/bin/SnapPy-x86_64.AppImage ~/bin/SnapPy
  ~/bin/SnapPy-x86_64.AppImage --install

The last command registers the SnapPy app with your desktop system and
starts it.  Next time, you can start SnapPy by using the desktop
search tool with GNOME or the main menu with KDE. You can pin also the icon
to your dash or task bar for easy access.  From a terminal window, you
can also start the app by typing ``SnapPy`` provided ``~/bin`` is in
your `$PATH <https://opensource.com/article/17/6/set-path-linux>`_.


Python Modules for macOS or Windows
-----------------------------------

If you write Python programs on macOS or Windows, you may wish to
install SnapPy as a Python module into your own copy of Python.  We
support Python 3.9 and up.  (On macOS, use a Python downloaded from
`Python.org <http://python.org>`_ and not the one provided by Apple.)
After installing Python, you may install a SnapPy module from your
Terminal application or Command Prompt with the command::

    python3 -m pip install --upgrade --user snappy snappy_15_knots

If you do not want the larger version of HTLinkExteriors that includes
the 15 crossing knots (it uses 110M of disk space), omit
``snappy_15_knots`` from the command.


Python Modules for Linux
------------------------

You can also use SnapPy with your Linux system's version of Python,
for example if you want to incorporate SnapPy in your own Python
scripts.  These instructions assume you have system administrator
(superuser) privileges to install software packages from your Linux
distribution.  (If you're not a superuser, use either the
`Linux app`_ or `conda`_.)

The first step is to install Python and other requirements; here's how
to do that on the most popular Linux distributions:

+ **Ubuntu/Debian/Mint/MX Linux/Elementary:** Tested on Ubuntu 24.04::

    sudo apt install python3-pip python3-tk

+ **Fedora**: Tested on Fedora 41::

    sudo dnf install python3-pip python3-tkinter

+ **Arch/Manjaro/EndeavourOS**: Install via::

    sudo pacman -Sy python-pip tk

+ **openSUSE**: For Tumbleweed::

    sudo zypper install python3-pip python3-tk

  For Leap, as of version 15.6 you need to ask for a recent version of
  Python or it will give you Python 3.6 which is too old for SnapPy::

    sudo zypper install python312-pip python312-tk

  You will need to replace ``python3`` by ``python3.12`` in subsequent
  steps.

+ **Red Hat Enterprise Linux/CentOS/Rocky Linux/AlmaLinux:**: These instructions
  are for version 8 or later; tested on AlmaLinux 8 and 9::

    sudo dnf install python3.11-pip python3.11-tkinter

  You will also need to replace ``python3`` by ``python3.11`` in subsequent
  steps.


Next, you need to install the SnapPy Python modules. For this, you can
either use a *virtual environment* (``python -m venv``) or do a *user
install* (``pip install --user``).  The former will work on any
version of Linux, whereas the latter is now strongly discouraged on
many systems (e.g. Ubuntu 24.04).  If you have not previously
installed SnapPy on this computer, we recommend using a virtual
environment, but suggest a user install if you are upgrading an
existing version of SnapPy that was installed in that manner.

Virtual environment
  Here is the `official tutorial
  <https://docs.python.org/3/tutorial/venv.html>`_ on using virtual
  environments in Python and an `in-depth article
  <https://realpython.com/python-virtual-environments-a-primer/>`_.  A
  recipe is::

    python3 -m venv snappy_venv
    # Switch to snappy_venv's Python
    source snappy_venv/bin/activate
    pip install snappy snappy_15_knots
    # Start the SnapPy app!
    SnapPy
    # Return to system Python
    deactivate

  If you always want to use the ``snappy_venv`` Python, adjust your
  `$PATH <https://opensource.com/article/17/6/set-path-linux>`_ to
  **start** with ``snappy_venv/bin``.

User install
  To do a user install with pip, try::

    # Note no "sudo" below!
    python3 -m pip install --upgrade --user snappy snappy_15_knots

  If you get a long error message that starts::

    error: externally-managed-environment

  you should probably use a virtual environment; however,
  you can force pip to do a user install via::

    # Note no "sudo" below!
    python3 -m pip install --upgrade --user --break-system-packages snappy snappy_15_knots

  Despite the scary name, provided you don't use ``sudo``, this will
  not actually modify the system packages, but rather install
  ``snappy`` into the subdirectory
  ``~/.local/share/python3.*/site-packages`` of your home directory.

  After a user install, you run the following command to start
  the app::

    ~/.local/bin/SnapPy

  So that you can start SnapPy with just the command ``SnapPy``, make
  sure ``~/.local/bin`` is in `in your path
  <https://opensource.com/article/17/6/set-path-linux>`_.


SageMath
--------

SnapPy has some special features when used within `SageMath
<http://sagemath.org>`_, the universal mathematics software based on
Python.  This section describes how to install SnapPy into your
existing copy of SageMath::

  sage -pip install --upgrade snappy snappy_15_knots

Alternatively, SageMath on `CoCalc <https://cocalc.com/>`_ (formerly
the SageMathCloud) also has SnapPy preinstalled, and the graphics
features even work via the `X11 interface
<http://blog.sagemath.com/cocalc/2018/11/05/x11.html>`_, see the
bottom of that page for more.

The graphical features may or may not work, depending on how Tkinter
was configured within Sage.  (There is no problem on macOS if you use
this `SageMath binary
<https://github.com/3-manifolds/Sage_macOS/releases>`_.)  If the
graphical features seem to "hang" when you try to start them, type
``%gui tk`` at the Sage prompt; please note that doing so may break
Sage's "attach" feature.


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
and many more; using this is sometimes the easiest way to get a
working setup for such multifaceted computations, especially on Windows.  For more, watch
`this demonstration <https://icerm.brown.edu/video_archive/?play=1992>`_.


Conda
-----

Conda can be used to install Python on all platforms and is a
particularly good choice to use SnapPy on the older Linux systems
often found on high-performance clusters.  Here is a recipe for
installing SnapPy into a new conda environment on macOS or Linux::

  source ~/miniforge3/bin/activate
  mamba create --name snappy_env python=3.12
  conda activate snappy_env
  pip install snappy
  python -m snappy.app


Source code
-----------

The complete source code for all platforms: |tarball|_

You can also browse our `source code repository
<https://github.com/3-manifolds/SnapPy>`_ or clone it using `git
<https://git-scm.com/>`_ via::

  git clone https://github.com/3-manifolds/SnapPy.git
