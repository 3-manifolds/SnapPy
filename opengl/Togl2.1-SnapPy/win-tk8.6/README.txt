This directory allows building the Togl extension for Windows using
Visual Studio 10.0.  This is the compiler used to build the python.org
distributions of Python 3.4.  It can be installed for free on a
Windows 7 system.  (The install will fail on Windows 10 because it is
impossible to install the required version 4.0 of the .NET Framework.)
Download links and instructions are here:
  https://wiki.python.org/moin/WindowsCompilers
Note that the link for the .NET Framework is not correct.  Instead use:
  https://www.microsoft.com/en-us/download/details.aspx?id=17851

The tcltk and tcltk64 subdirectories contain the needed headers and
libraries from the python.org distribution of Python 2.7.9.

To build Togl 2.1:
  * Open a "command prompt". (Press <Windows-X> or run "cmd" in the Msys shell.)
  * Change to this directory. (Running cmd in Msys preserves the cwd.)
  * To build for 32-bit:
    >make-togl.bat /x86
  * The Togl 2.0 installation will be in win32VC-tk8.6
  * To build for 64-bit:
    >make-togl.bat /x64
  * The Togl 2.0 installation will be in win32VC-x86_64-tk8.6
  * To remove intermediate files:
    >make-togl.bat /x86 clean
    >make-togl.bat /x64 clean
  * To remove all generated files:
    >make-togl.bat /x86 distclean
    >make-togl.bat /x64 distclean
  
