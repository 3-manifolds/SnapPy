Warmup
--------------

1. Upload plink, spherogram, and snappy tarballs to testpypi. 

2. Test source tarballs on Linux build boxes via::

     cd SnapPy/dev/release_tools
     py27 test_pypi -p -t snappy
     py26 test_pypi -p -t snappy

   Fix all Python 2.6 incompatible syntax that has gotten into the
   packages.  Upload new tarballs to testpypi, and retest.  

3. Build Mac disk image on old Mac Pro (10.6).  Test on new Mac Pro (10.9)
   and Mac Book (10.10).

4. Build Window exe on Win7, test on Windows 8.1.

   


   

     sage -pip install --no-cache-dir -i https://testpypi.python.org/simple snappy --upgrade --user


     
