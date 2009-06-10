.. Installing SnapPy

Installing and running SnapPy
======================================================

Here are detailed instructions on how to get SnapPy working on a
variety of platforms.

Mac OS X
---------------

Simply download `SnapPy.dmg <http://math.uic.edu/~t3m/SnapPy/SnapPy.dmg>`_
and copy SnapPy.app to the Applications folder.  Double-click to start
it, just like any other application.


Linux
--------------------

We will eventually have some easy-to-install binaries here, but for
now, follow the instructions for `generic Unix`_ in the next section.

Generic Unix
----------------------------------------------------------

If you use a Unix other that OS X or Linux, or if the prebuilt
packages don't work for you, you'll need to build SnapPy from source.
Here are some detailed instructions.

Things you'll need:

- Python with Tkinter: You'll need to have `Python
  <http://python.org>`_ (at least version 2.5) and `Tk <http://tcl.tk>`_
  (at least version 8.4) with `Tkinter <http://wiki.python.org/moin/TkInter>`_ to
  connect them, including the header files.  For instance, on Debian or
  Ubuntu, install the packages "python-tk" and "python-dev".  In
  addition, you'll need

  - `Setuptools <http://peak.telecommunity.com/DevCenter/setuptools>`_,
    which is typically packaged as "python-setuptools" or can be installed via::

      curl -O http://peak.telecommunity.com/dist/ez_setup.py
      sudo python ez_setup.py  

  Test that Python is in order by installing PLink from source::

    curl -O http://math.uic.edu/~t3m/plink/plink.tar.gz
    tar xfz plink.tar.gz; cd plink
    sudo python setup.py install 
    cd /tmp; python -m plink.app   # Should start the link editor!

- Support for OpenGL (3D graphics): This is built in on OS X, but on other Unixes,
  you'll need to install `MESA <http://www.mesa3d.org/>`_ and `FreeGLUT
  <http://freeglut.sf.net>`_.  For instance on Debian 
  and Ubuntu, install "freeglut3-dev".  

- `Cython <http://cython.org>`_, which you can install via::

    sudo python -m easy_install cython

Now download the `source code`_ listed below, for instance::

    curl -O http://math.uic.edu/~t3m/SnapPy/SnapPy.tar.gz
    tar xfz SnapPy.tar.gz; cd SnapPy

There are two more dependencies that need to be dealt with:

- `PARI <http://pari.math.u-bordeaux.fr/>`_:  Inside the SnapPy directory do::

    bash build_pari.sh   # Downloads and builds the PARI library

  
- `Togl <http://togl.sf.net>`_: a 3d widget for Tk. For OS X and
  Linux, there are pre-built binaries of this in the snappy
  subdirectory, e.g. snappy/linux2-tk8.4.  First, test if those work
  via::

    sudo python -m easy_install PyOpenGL
    python snappy/polyviewer.py     

  A Dirichlet domain should appear, which you can spin around etc. If
  this doesn't work, you'll need to edit or follow "build_togl.sh" to
  build Togl directly, and create the appropriate subdirectory of
  snappy.

  
Finally, compile and install the SnapPy module (which will install
certain other dependencies) and test::

  sudo python setup.py install
  cd /tmp; python -m snappy.app

You may get a message about creating a ".ipython" directory; this is
normal, just hit return to continue.  There should also now be a
command "snappy" which does the same thing as "python -m snappy.app".

Sage
----

SnapPy has some special features when used within `Sage
<http://sagemath.org>`_, the universal mathematics software based on
Python. Installation is easy::

 curl -O http://math.uic.edu/~t3m/SnapPy/SnapPy.tar.gz
 tar xfz SnapPy.tar.gz; cd SnapPy
 sage -python setup.py install

The graphics features may or may not work, depending on how Tkinter
was configured within Sage, but everything else will work fine.

Windows
-------------------

Not yet available, though we plan on this in future. If you're familiar
with `py2exe <http://py2exe.org>`_ and MVC feel free to get this
working for us.


Source code
-----------------------------------

The complete source code for all platforms: `SnapPy.tar.gz <http://math.uic.edu/~t3m/SnapPy/SnapPy.tar.gz>`_   

You can also get it straight from the `Mercurial
<www.selenic.com/mercurial>`_ repository::

  hg clone static-http://math.uic.edu/~t3m/hg/SnapPy





