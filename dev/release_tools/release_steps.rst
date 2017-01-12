Overview
========

Key tools:

1. macOS binaries are built on a OS X 10.6 VM running in VirtualBox on
   thurston.math.illinois.edu.

2. Linux binaries are build via a derivative of the manylinux1 Docker
   image; for details, see "docker/README.rst".

3. Windows wheels are built automatically via AppVeyor, and can be
   accessed via the "Artifacts" tab on the job page.

4. The Windows app is build on a Win10 64bit VM running in VMWare
   Fusion on thurston.math.illinois.edu.

5. The script "test_pypi.py" is a key tool. It creates a virtual
   environment for testing a package posted on (test)pypi.python.org.

6. Linux testing done on a Ubuntu 16.04 (64-bit) VM running in VMWare
   Fusion on thurston.math.illinois.edu.
   

Warmup
======

1. Upload plink, spherogram, and snappy tarballs to testpypi. Be sure
   to use "rc*" version names since it won't let you use the same name
   twice, even on this test server::

     cd blah
     python setup.py sdist
     twine upload -r test dist/blah.tar.gz

   Further details can be found in "pypi.rst".

2. Fire up the Linux testing VM and do::

     cd SnapPy/dev/release_tools
     hg pull -u
     py27 test_pypi.py -p -t snappy
     py35 test_pypi.py -p -t snappy

3. Build Mac disk image and wheel 10.6 VM.  Test on that machine and
   some newer one as well by via starting the app and typing::

     import snappy.test

4. Build Window exe on Win7, test on that machine and Win10 via
   installing the app and typing::

     import snappy.test

5. Do doctests in Sage.


6. **Future** Contact beta testers, providing incantations::

     python -m pip install --pre --extra-index-url https://testpypi.python.org/simple --upgrade --no-deps plink spherogram snappy

   or for those who use sage::

       sage -pip install --pre --extra-index-url https://testpypi.python.org/simple --upgrade --no-deps --no-binary :all: plink spherogram snappy

7. Fix issues pointed out by beta testers.  Lather, rinse, repeat.

   Possible beta testers: Ken Baker, Craig Hodgson, Dave Futer, Saul
   Schleimer, Mark Bell, and Ilya Kofman.

Actual release
----------------------

1. Remove "rc" from "version.py", **rebuild snappy docs, check
   results**, and *commit and push*.  Then rebuild tarballs.

2. Build the eggs and wheels for all platforms, putting the results in
   some directory.  On linux, use the "update_SnapPy_linux.py" script.

3. Rebuild full apps for OS X and Windows.  Test one last time.

4. Use twine to upload everything to PyPI.

5. Upload OS X and Windows apps to Bitbucket.  Record downloads. 

6. Update documentation on web.

7. Update "current.txt".


Final testing
=============

On various systems, do::

  cd SnapPy/dev/release_tools
  py27 test_pypi.py -p/-e snappy

including SageMathCloud (but there's an issue with sys.path in the
latter).  When testing the stand-alone apps, do "import snappy.test".

Then tag the releases in Mercurial::

  hg tag 1.4_as_released; hg push




Announce to the world
=====================

1. LDT blog

2. Facebook

3. Google+

4. Mailing list

5. William Stein 


Application Download Counts
===========================

a. Version 2.3.*: 796 Mac, 955 Windows.
b. Version 2.4.0:  16 Mac, 18 Windows.
c. Version 2.4.1:  352 Mac, 935 Windows.


