@echo off

@rem Process command line args
if "%1" =="-arch" (
   set arch=%2
   shift & shift
   )
if "%arch%"=="" set arch=x86
if "%1"=="" (set target=Togl) else (set target=%1)

@rem This is the standard install directory for VC for Python27
@rem If you move it, change this.
set VCFORPYTHONLONG=C:\Program Files (x86)\Common Files\Microsoft\Visual C++ for Python\9.0

@rem Build the short path to the VC directory
set VCFORPYTHON=
for %%P in ("%VCFORPYTHONLONG%") do call :extendpath %%~sP
set VCFORPYTHONSDK=%VCFORPYTHON%\WinSDK

@rem The 32-bit compiler
set cc32=%VCFORPYTHON%\VC\bin\cl.exe

@rem The Resource Compiler
set RC=%VCFORPYTHONSDK%\Bin\RC

@rem To avoid surprises we set a minimal path.
set PATH=C:\WINDOWS\system32;C:\WINDOWS;C:\WINDOWS\System32\Wbem;

@rem Switch, depending on architecture.
if /i "%arch%"=="x86" goto :x86
if /i "%arch%"=="x64" goto :x64
goto :usage

:x86
echo Setting the environment for a 32 bit x86 build.
set CCDIR=%VCFORPYTHON%\VC\bin
if not exist "%CCDIR%"\cl.exe goto missing
@rem The vcbuild tool is a 64 bit application. So we add bin\amd64 to the path. 
set PATH=%VCFORPYTHON%\VC\bin;%VCFORPYTHONSDK%\Bin;%VCFORPYTHON%\VC\bin\amd64;%PATH%
set INCLUDE=%VCFORPYTHON%\VC\include;%VCFORPYTHONSDK%\Include;
set LIB=%VCFORPYTHON%\VC\lib;%VCFORPYTHONSDK%\Lib;
set LIBPATH=%VCFORPYTHON%\VC\lib;%VCFORPYTHONSDK%\Lib;
set CC=%VCFORPYTHON%\VC\bin\cl.exe
set TCLDIR=tcltk
set MACHINE=IX86
set TARGET_CPU=x86
set CPUTAG=
goto :run
)

:x64
echo Setting the environment for a 64 bit build.
set CCDIR=%VCFORPYTHON%\VC\bin\amd64
if not exist "%CCDIR%"\cl.exe goto missing
set PATH=%VCFORPYTHON%\VC\bin\amd64;%VCFORPYTHONSDK%\Bin\x64;%PATH%
set INCLUDE=%VCFORPYTHON%\VC\include;%VCFORPYTHONSDK%\Include;
set LIB=%VCFORPYTHON%\VC\lib\amd64;%VCFORPYTHONSDK%\Lib\x64;
set LIBPATH=%VCFORPYTHON%\VC\lib\amd64;%VCFORPYTHONSDK%\Lib\x64;
set CC=%VCFORPYTHON%\VC\bin\amd64\cl.exe
set TCLDIR=tcltk64
set MACHINE=X64
set TARGET_CPU=x86_64
set CPUTAG=-x86_64
goto :run
)

:missing
echo The %1 compiler was not found in %CCDIR%
echo Check that the install path is correct in setup-env.bat.
goto :eof

:usage
echo got %arch% %target% %arg%
echo usage: %~0 [-arch x86^|x64] [Togl^|clean^|veryclean^|distclean]
echo Both args are optional and can come in either order.
goto :eof

:run
echo Current Environment:
path
echo ==
echo INCLUDE=%INCLUDE%
echo ==
echo LIB=%LIB%
echo ==
echo VCFORPYTHON=%VCFORPYTHON%
echo ==
echo CC=%CC%
nmake -f makefile.vc %target%
if /i "%target%"=="Togl" (
  nmake -f makefile.vc install TARGET_CPU=%TARGET_CPU% INSTALLDIR=win32VC%CPUTAG%-tk8.5
)

:extendpath
set VCFORPYTHON=%VCFORPYTHON%%1
goto :eof
