Overview
========

Key tools:

1. Wheels for all three platforms are built using cibuildwheel as part
   of an automatic GitHub action after each commit.

2. The sdist tarballs are built by a GitHub action of the
   "snappy_release" project.

3. The Windows and macOS apps is are built by a GitHub actions of
   "snappy_release".

4. The script "test_pypi.py" is a key tool. It creates a virtual
   environment for testing a package posted on (test)pypi.python.org.

5. Testing done on virtual machines running in VMWare Fusion on
   dehn.math.illinois.edu.

6. Nathan stores current and old versions in "~/Dropbox/pypi/".


Warmup
======

0. Bump versions numbers for each subproject to whatever you want them
   to be for final release, plus some "rcN" suffix.  Commit and push.

1. Trigger the GitHub actions for "snappy_release".  Download the
   artifacts for the "sdist" action and upload to testpypi::

      twine upload -r test dist/snappy/*.tar.gz

   Further details can be found in "pypi.rst".

2. Fire up a Linux test VM and do::

     py36 test_pypi.py -p -t snappy

   Repeat for macOS and Windows.  This will flush out any new issues
   with the sdist tarballs.

3. Download Windows installer from "snappy_release" actions and test.

4. Do doctests in the most recent beta Sage if needed (handled by the
   CI for snappy *if* the Sage docker images are up to date).

5. Upload GH built wheels to testpypi and test.  Include older
   versions of macOS.


Actual release
==============

1. Bump version numbers to final, trigger "snappy_release" rebuild.

2. Build final macOS app.

3. Use twine to upload everything to PyPI.

4. Create releases on GitHub and upload macOS and Windows apps.

5. Update documentation on web by copying docs into t3m_web/SnapPy.

6. **Create** "current.txt" in t3m_web/SnapPy.


Final testing
=============

On various systems, do::

  cd SnapPy/dev/release_tools
  py27 test_pypi.py -p/-e snappy

including CoCalc (but there's an issue with sys.path in the
latter).  When testing the stand-alone apps, do "import snappy.test;
snappy.test.runtests()"

A super-fast way to check the non-graphical stuff on Linux using
Docker and the official Python images is::

  docker run -it python:3.4-stretch /bin/bash
  pip install snappy; python -m snappy.test



Announce to the world
=====================

1. Facebook

2. Mailing list

3. CoCalc ticket


Application Download Counts
===========================

a. Version 2.3.*:  796 Mac,  955 Windows.
b. Version 2.4.*:  471 Mac, 1048 Windows.
c. Version 2.5.*:  433 Mac,  729 Windows.
d. Version 2.6.0:  383 Mac,  699 Windows (28% Python 3).
e. Version 2.6.1: 1018 Mac, 1129 Windows (15% Python 3).
f. Version 2.7:    986 Mac,  636 Windows (89% Python 3).
g. Version 2.8:   295* Mac,  573 Windows (85% Python 3).

   Mac app was silently replaced a couple times to deal with Big Sur
   and Tk issues.  So missing downloads from Oct 11 through Feb 21; 295
   is estimate based on 175 known downloads.

h. Version 3.0(.1): 292 Mac, 434 Windows.
i. Version 3.0.2:    68 Mac,  70 Windows.
j. Version 3.0.3:   350 Mac  645 Windows. (Through Oct 11)

Average downloads for 2015-3-22 through 2017-10-26.

Mac (app on bitbucket): 53/month
Windows (app on bitbucket): 86/month
PyPI (eggs/wheels/tarballs): 161/month

PyPI statistics are based on:

http://www.pypi-stats.com/package/?q=snappy

for just the period of 2017.


Average downloads for 2017-11-27 through 2021-10-4 (46 months).

Mac (app on bitbucket/github): 63/month
Windows (app on bitbucket/github): 75/month
Docker (kitchen sink): 44/month
PyPI (eggs/wheels/tarballs): ??/month

The PyPI numbers are so large as to be unbelievable, more than 10,000
a month. Nearly all of these are for Linux, with "just" 600 a month
for Windows and Mac combined.




Getting download stats from GitHub:

https://api.github.com/repos/3-manifolds/SnapPy/releases
