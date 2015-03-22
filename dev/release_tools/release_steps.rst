Warmup
--------------

1. Upload plink, spherogram, and snappy tarballs to testpypi. Be sure
   to use "rc" version names since it won't let you use the same name
   twice, even on this test server.

2. Test source tarballs on Linux build boxes via::

     cd SnapPy/dev/release_tools
     py27 test_pypi.py -p -t snappy
     py26 test_pypi.py -p -t snappy

   Fix all Python 2.6 incompatible syntax that has gotten into the
   packages.  Upload new tarballs to testpypi, and retest.  

3. Build Mac disk image and wheel on old Mac Pro (10.6).  Test on new Mac Pro (10.9)
   and Mac Book (10.10).

4. Build Window exe on Win7, test on Windows 8.1.

Actual release
----------------------

1. Remove "rc" from "version.py", **rebuild snappy docs**, and rebuild
   tarballs.

2. Build the eggs and wheels for all platforms, putting the results in
   some directory.

3. Rebuild full apps for OS X and Windows.  Test one last time.

4. Use twine to upload everything to PyPI.

5. Upload OS X and Windows apps to Bitbucket.

6. Update documentation on web.

7. Update "current.txt".


Final testing:
---------------------------

On various systems, do::

  cd SnapPy/dev/release_tools
  py27 test_pypi.py -p/-e snappy

including SageMathCloud (but there's an issue with sys.path in the
latter).

Then tag the releases in Mercurial::

  hg tag 1.4_as_released; hg push


Announce to the world:
---------------------------

1. LDT blog

2. Facebook

3. Google+

4. Mailing list

5. William Stein 
