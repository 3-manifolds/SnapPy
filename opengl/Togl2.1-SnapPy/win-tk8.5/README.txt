This directory allows building the Togl extension for Windows using
the "Visual C++ for Python" package, which can be downloaded for free
from Microsoft:

http://www.microsoft.com/en-us/download/details.aspx?id=44266

The tcltk and tcltk64 subdirectories contain the needed headers and
libraries from the python.org distribution of Python 2.7.9.

To build Togl 2.1:
  * Open a "command prompt". (Press <Windows-X> or run "cmd" in the Msys shell.)
  * Change to this directory. (Running cmd in Msys preserves the cwd.)
  * run:
    >make-togl.bat
  * The Togl 2.0 installation will be in win32VC-tk8.5
  * To build for 64-bit:
    >make-togl.bat -arch x64
  * To remove intermediate files:
    >make-togl.bat clean
  * To remove all generated files:
    >make-togl.bat distclean
  
