@echo off
if "%1" == "" goto x86
if not "%2" == "" goto usage

if /i %1 == x86       goto x86
if /i %1 == amd64     goto amd64
if /i %1 == x64       goto amd64
if /i %1 == ia64      goto ia64
if /i %1 == x86_amd64 goto x86_amd64
if /i %1 == x86_ia64  goto x86_ia64
goto usage

:x86
call "C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv.Cmd" /Release /x86
goto :eof

:amd64
call "C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv.Cmd" /Release /x64
goto :eof

:ia64
if not exist "%~dp0bin\ia64\vcvars64.bat" goto missing
call "%~dp0bin\ia64\vcvars64.bat"
goto :eof

:x86_amd64
if not exist "%~dp0bin\x86_amd64\vcvarsx86_amd64.bat" goto missing
call "%~dp0bin\x86_amd64\vcvarsx86_amd64.bat"
goto :eof

:x86_ia64
if not exist "%~dp0bin\x86_ia64\vcvarsx86_ia64.bat" goto missing
call "%~dp0bin\x86_ia64\vcvarsx86_ia64.bat"
goto :eof

:usage
echo Error in script usage. The correct usage is:
echo     %0 [option]
echo where [option] is: x86 ^| ia64 ^| amd64 ^| x86_amd64 ^| x86_ia64
echo:
echo For example:
echo     %0 x86_ia64
goto :eof

:missing
echo The specified configuration type is missing.  The tools for the
echo configuration might not be installed.
goto :eof
