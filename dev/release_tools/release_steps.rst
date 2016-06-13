Warmup
--------------

1. Upload plink, spherogram, and snappy tarballs to testpypi. Be sure
   to use "rc" version names since it won't let you use the same name
   twice, even on this test server::

     twine upload -r test dist/blah.tar.gz

2. Test source tarballs on Linux build boxes via::

     cd SnapPy/dev/release_tools
     hg pull -u; hg branch; hg status
     py27 test_pypi.py -p -t snappy
     py26 test_pypi.py -p -t snappy

   Fix all Python 2.6 incompatible syntax that has gotten into the
   packages.  Upload new tarballs to testpypi, and retest.  

3. Build Mac disk image and wheel Mac Pro VM running Snow Leopard
   10.6.  Test that machine, the new Mac Pro (10.9), and the Mac Book
   (10.11) via starting the app and typeing::

     import snappy.test

4. Build Window exe on Win7, test on that machine and Win10 via
   installing the app and typing::

     import snappy.test

5. Do doctests in Sage.


6. Contact beta testers, providing incantations::

     python -m pip install --pre --extra-index-url https://testpypi.python.org/simple --upgrade --no-deps plink spherogram snappy

     or for those who use sage::

       sage -pip install --pre --extra-index-url https://testpypi.python.org/simple --upgrade --no-deps --no-use-wheel plink spherogram snappy

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


Final testing:
---------------------------

On various systems, do::

  cd SnapPy/dev/release_tools
  py27 test_pypi.py -p/-e snappy

including SageMathCloud (but there's an issue with sys.path in the
latter).  When testing the stand-alone apps, do "import snappy.test".

Then tag the releases in Mercurial::

  hg tag 1.4_as_released; hg push




Announce to the world:
---------------------------

1. LDT blog

2. Facebook

3. Google+

4. Mailing list

5. William Stein 


Application Download Counts:
-------------------------------------

a. Version 2.3.*: 796 Mac, 955 Windows.
b. Version 2.4.1:  16 Mac, 18 Windows.


