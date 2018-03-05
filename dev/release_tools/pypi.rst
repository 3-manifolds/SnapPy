Notes on using PyPI
===================

To start, make sure you have the latest versions of some key packages::

  python -m pip install --upgrade pip setuptools twine virtualenv

With Python 2.7.9 or Python 3.4, you can install pip via::

  python -m ensurepip

if it is missing.


Building the source package
===========================

The first step is to build the tarball and test on a local machine::

  python setup.py clean sdist

which creates "dist/package-version.tar.gz".  Now create a virtualenv
in the "pt" directory for testing::

  rm -rf pt; virtualenv pt
  pt/bin/pip install --upgrade setuptools

and install from the source distribution::

  pt/bin/pip install dist/package-version.tar.gz
  pt/bin/python -c 'import package; print package.__version__'
  pt/bin/python -m package.test -v 

Hints: May need to create a "MANIFEST.in" file to get everything put
in the "sdist" bundle correctly.  

  
Uploading a source-only package toPyPI
======================================

Initially, it is best to work with the sandbox version of PyPI::

  https://testpypi.python.org/pypi

Here are the steps to create a source-only package on TestPyPI.

a. Go to::

     https://testpypi.python.org/pypi?%3Aaction=submit_form

   Fill out only the project name and version number "0.0" and then
   click the "Add information" button at the bottom.

b. Assuming one's ".pypirc" is correctly configured, just do::

   twine upload -r test dist/package-version.tar.gz

c. Now test it::

   rm -rf pt; virtualenv pt
   pt/bin/pip install --no-cache-dir -i https://testpypi.python.org/simple package
   pt/bin/python -m package.test -v

Now to do the real thing.

i. Go to https://pypi.python.org/pypi?%3Aaction=submit_form and repeat
   step (a) above.

ii. Upload, then test::

      twine upload dist/package-version.tar.gz
      


   
