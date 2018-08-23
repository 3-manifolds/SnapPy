Overview
========

Key tools:

1. macOS binaries are built on a OS X 10.12 VM running in VMWare Fusion on
   thurston.math.illinois.edu.

2. Linux binaries are build via a derivative of the manylinux1 Docker
   image; for details, see "docker/README.rst".

3. Windows wheels are built automatically via AppVeyor, and can be
   accessed via the "Artifacts" tab on the job page.

4. The Windows app is built via AppVeyor and the "snappy_release" project.

5. The script "test_pypi.py" is a key tool. It creates a virtual
   environment for testing a package posted on (test)pypi.python.org.

6. Linux testing done on a Ubuntu 18.04 (64-bit) VM running in VMWare
   Fusion on thurston.math.illinois.edu.

7. Nathan stores current and old versions in "~/Dropbox/pypi/".
   

Warmup
======

0.  Bump versions numbers to whatever you want them to be for the
    final release.  Commit and push.  

1.  For each of plink, sphereogram, and snappy, run "python setup.py
    release" on one platform.  This will generate, among other things,
    an "sdist" tarball.  Since even testpypi will not allow you to use
    the same name twice, add a release candidate tag and then upload
    to testpypi::

      rctag.py -r1 dist/*.tar.gz
      twine upload -r test dist/*.tar.gz

   Further details can be found in "pypi.rst".

2. Fire up the Linux testing VM and do::

     cd SnapPy/dev/release_tools
     hg pull -u
     py27 test_pypi.py -p -t snappy
     py35 test_pypi.py -p -t snappy

3. Build Mac disk image and wheel on 10.12 VM.  Test on that machine and
   some older one as well by via starting the app and typing::

     import snappy.test
     snappy.test.runtests()

4. Build Window exe on Win10, test on that machine and Win7 via
   installing the app and typing::

     import snappy.test
     snappy.test.runtests()

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

1. Build the eggs and wheels for all platforms, putting the results in
   some directory.  

2. Rebuild full apps for OS X and Windows.  Test one last time.

3. Use twine to upload everything to PyPI.

4. Upload OS X and Windows apps to Bitbucket.  Record downloads. 

5. Update documentation on web.

6. Update "current.txt".


Final testing
=============

On various systems, do::

  cd SnapPy/dev/release_tools
  py27 test_pypi.py -p/-e snappy

including SageMathCloud (but there's an issue with sys.path in the
latter).  When testing the stand-alone apps, do "import snappy.test".

Then tag the releases in Mercurial::

  hg tag 1.4_as_released; hg push

A super-fast way to check the non-graphical stuff on Linux using
Docker and the official Python images is::

  docker run -it python:3.4-stretch /bin/bash
  pip install snappy; python -m snappy.test



Announce to the world
=====================

1. LDT blog

2. Facebook

3. Google+

4. Mailing list

5. William Stein 


Application Download Counts
===========================

a. Version 2.3.*: 796 Mac,  955 Windows.
b. Version 2.4.*: 471 Mac, 1048 Windows.
c. Version 2.5.*: 433 Mac, 729 Windows.
d. Version 2.6.0: 383 Mac, 699 Windows (28% Python 3)
e. Version 2.6.1: 


Average downloads for 2015-3-22 through 2017-10-26.

Mac (app on bitbucket): 53/month
Windows (app on bitbucket): 86/month
PyPI (eggs/wheels/tarballs): 161/month

PyPI statistics are based on:

http://www.pypi-stats.com/package/?q=snappy

for just the period of 2017.  

