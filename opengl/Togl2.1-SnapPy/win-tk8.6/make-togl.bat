@ECHO OFF
:: --------------------------------------------------------------------------------------------
:: File : SetEnv.cmd
::
:: Abstract: This batch file sets the appropriate environment variables for the Windows SDK
::           build environment with respect to OS and platform type.
::
:: Usage : Setenv [/Debug | /Release][/x86 | /x64 | /ia64 ][/vista | /xp | /2003 | /2008 | /win7][-h | /?]
::
::                /Debug   - Create a Debug configuration build environment
::                /Release - Create a Release configuration build environment
::                /x86     - Create 32-bit x86 applications
::                /x64     - Create 64-bit x64 applications
::                /ia64    - Create 64-bit ia64 applications
::                /vista   - Windows Vista applications
::                /xp      - Create Windows XP SP2 applications
::                /2003    - Create Windows Server 2003 applications
::                /2008    - Create Windows Server 2008 or Vista SP1 applications
::                /win7    - Create Windows 7 applications
::
:: Note: This file makes use of delayed expansion to get around the limitations
:: of the current expansion which happens when a line of text is read , not when it is executed.
:: For more information about delayed expansion type: SET /? from the command line
:: --------------------------------------------------------------------------------------------


:: --------------------------------------------------------------------------------------------
:: Set the default value for target and current CPU based on processor architecture.
:: --------------------------------------------------------------------------------------------
::SET "CURRENT_CPU=x86"
::SET "TARGET_CPU=x86"
SET "Configuration=Release"
::SET "TARGET_PLATFORM=WIN7"

IF "x%TARGET_CPU%x"=="xx" (
IF /I "%PROCESSOR_ARCHITECTURE%"=="x86" SET "TARGET_CPU=x86" & SET "CURRENT_CPU=x86"
IF /I "%PROCESSOR_ARCHITEW6432%"=="x86" SET "TARGET_CPU=x86" & SET "CURRENT_CPU=x86"
IF /I "%PROCESSOR_ARCHITECTURE%"=="AMD64" SET "TARGET_CPU=x64" & SET "CURRENT_CPU=x64"
IF /I "%PROCESSOR_ARCHITEW6432%"=="AMD64" SET "TARGET_CPU=x64" & SET "CURRENT_CPU=x64"
IF /I "%PROCESSOR_ARCHITECTURE%"=="x64"   SET "TARGET_CPU=x64" & SET "CURRENT_CPU=x64"
IF /I "%PROCESSOR_ARCHITECTURE%"=="IA64"  SET "TARGET_CPU=IA64" & SET "CURRENT_CPU=IA64"
IF /I "%PROCESSOR_ARCHITEW6432%"=="IA64"  SET "TARGET_CPU=IA64" & SET "CURRENT_CPU=IA64"
GOTO Parse_Args
)
:: --------------------------------------------------------------------------------------------
:: Parse command line argument values.
:: Note: For ambiguous arguments the last one wins (ex: /debug /release)
:: --------------------------------------------------------------------------------------------
:Parse_Args
IF /I "%~1"=="/debug"      SET "Configuration=Debug"   & SHIFT & GOTO Parse_Args
IF /I "%~1"=="/release"    SET "Configuration=Release" & SHIFT & GOTO Parse_Args
IF /I "%~1"=="/x86"        SET "TARGET_CPU=x86"           & SHIFT & GOTO Parse_Args
IF /I "%~1"=="/x64"        SET "TARGET_CPU=x64"           & SHIFT & GOTO Parse_Args
IF /I "%~1"=="/ia64"       SET "TARGET_CPU=IA64"          & SHIFT & GOTO Parse_Args
IF /I "%~1"=="/vista"      SET "TARGET_PLATFORM=LH"       & SHIFT & GOTO Parse_Args
IF /I "%~1"=="/xp"         SET "TARGET_PLATFORM=XP"       & SHIFT & GOTO Parse_Args
IF /I "%~1"=="/2003"       SET "TARGET_PLATFORM=SRV"      & SHIFT & GOTO Parse_Args
IF /I "%~1"=="/2008"       SET "TARGET_PLATFORM=LHS"      & SHIFT & GOTO Parse_Args
IF /I "%~1"=="/win7"       SET "TARGET_PLATFORM=WIN7"     & SHIFT & GOTO Parse_Args
IF /I "%~1"=="-h"          GOTO Error_Usage
IF /I "%~1"=="/?"          GOTO Error_Usage
IF    "x%~1"=="x"          GOTO Done_Args

ECHO Unknown command-line switch: %~1
GOTO Error_Usage

:Done_Args

:: --------------------------------------------------------------------------------------------
:: Prevent path duplication if setenv is run multiple times within a single command session
:: --------------------------------------------------------------------------------------------
IF "x%ORIGINALPATH%x"=="xx" (
  SET "ORIGINALPATH=%PATH%"
) ELSE (
  SET "PATH=%ORIGINALPATH%"
)

:: --------------------------------------------------------------------------------------------
:: Default the build configuration to DEBUG if one is not specified.
:: Set the command prompt text color based on the build configuration.
:: --------------------------------------------------------------------------------------------
IF "x%Configuration%"=="x" SET Configuration=Debug
IF "%Configuration%"=="Debug" (
    SET DEBUG=1
    SET DEBUGMSG=Debug
        COLOR 03
) ELSE IF "%Configuration%"=="Release" (
    SET DEBUG=0
    SET DEBUGMSG=Release
        COLOR 02
) ELSE GOTO Error_Usage

:: --------------------------------------------------------------------------------------------
:: Default to LHS if no target type specified and configure for appropriate target
:: --------------------------------------------------------------------------------------------
IF "x%TARGET_PLATFORM%"=="x" (
FOR /F "tokens=1,2,3 delims=;." %%i IN ('Cmd /c Ver') DO (
    IF "%%i"=="Microsoft Windows XP [Version 5" (
      SET TARGET_PLATFORM=XP
    ) ELSE IF "%%i"=="Microsoft Windows [Version 5" (
      SET TARGET_PLATFORM=SRV
    ) ELSE IF "%%i"=="Microsoft Windows [Version 6" (
      IF "%%k" == "6000]" (
        SET TARGET_PLATFORM=LH
      ) ELSE IF "%%k" == "6001]" (
        SET TARGET_PLATFORM=LHS
      ) ELSE IF "%%k" == "7600]" (
        SET TARGET_PLATFORM=WIN7
    ) ELSE (
        SET TARGET_PLATFORM=WIN7
    )
    ) ELSE (
      SET TARGET_PLATFORM=WIN7
    )
  )
)

:: --------------------------------------------------------------------------------------------
:: Set the PlatformToolset that is used by MSBuild 4.0
:: To use the Visual Studio 2010 tools/headers/libraries, "set PlatformToolset=100"
:: --------------------------------------------------------------------------------------------

SET PlatformToolset=Windows7.1SDK
SET ToolsVersion=4.0
SET WindowsSDKVersionOverride=v7.1

:: --------------------------------------------------------------------------------------------
:: Determine which registry keys to look at based on architecture type.  Set default values for
:: VC and VS root, which would be used if one or both the corresponding registry keys are not
:: found.
:: --------------------------------------------------------------------------------------------
SET RegKeyPath=HKEY_LOCAL_MACHINE\SOFTWARE\Microsoft\VisualStudio\SxS\VC7
SET VSRegKeyPath=HKEY_LOCAL_MACHINE\SOFTWARE\Microsoft\VisualStudio\SxS\VS7
SET WinSDKRegKeyPath=HKEY_LOCAL_MACHINE\SOFTWARE\Microsoft\Microsoft SDKs\Windows\v7.1
:: find platform for default paths
IF "%CURRENT_CPU%"=="x86" (
  SET Path32=%ProgramFiles%
) ELSE (
  SET Path32=%ProgramFiles(x86)%
)
:: set VS default paths
SET "VCINSTALLDIR=%Path32%\Microsoft Visual Studio 10.0\VC\"
SET "VSINSTALLDIR=%Path32%\Microsoft Visual Studio 10.0\"

:: clear the temp variable
SET Path32=

:: set WinSDK default dir
SET WindowsSDKDir=%ProgramFiles%\Microsoft SDKs\Windows\v7.1\

IF EXIST "%WinDir%\Microsoft.NET\Framework\msbuild.exe" SET "FrameworkDir32=%WinDir%\Microsoft.NET\Framework\"
IF EXIST "%WinDir%\Microsoft.NET\Framework64\msbuild.exe" SET "FrameworkDir64=%WinDir%\Microsoft.NET\Framework64"

IF EXIST "%WinDir%\Microsoft.NET\Framework\v3.5\MSBuild.exe" SET "Framework35Version=v3.5"

:: Set the WindowsSdkDir path
FOR /F "tokens=2* delims=	 " %%A IN ('REG QUERY "%WinSDKRegKeyPath%" /v InstallationFolder') DO SET WindowsSDKDir=%%B
SET "sdkdir=%WindowsSDKDir%"

:: Set the framework paths and versions
FOR /F "tokens=2* delims=	 " %%A IN ('REG QUERY "%RegKeyPath%" /v FrameworkDir32') DO SET FrameworkDir32=%%B
FOR /F "tokens=2* delims=	 " %%A IN ('REG QUERY "%RegKeyPath%" /v FrameworkDir64') DO SET FrameworkDir64=%%B
FOR /F "tokens=2* delims=	 " %%A IN ('REG QUERY "%RegKeyPath%" /v FrameworkVer32') DO SET FrameworkVersion32=%%B
FOR /F "tokens=2* delims=	 " %%A IN ('REG QUERY "%RegKeyPath%" /v FrameworkVer64') DO SET FrameworkVersion64=%%B

SET FrameworkDir32=%WinDir%\Microsoft.NET\Framework\
SET "FrameworkDir=%FrameworkDir32%"
SET "FrameworkVersion=%FrameworkVersion32%"

IF NOT "%CURRENT_CPU%"=="x86"  GOTO SetRegPathFor64Bit

GOTO Done_SetRegPath

:SetRegPathFor64Bit
SET RegKeyPath=HKEY_LOCAL_MACHINE\SOFTWARE\Wow6432Node\Microsoft\VisualStudio\SxS\VC7
SET VSRegKeyPath=HKEY_LOCAL_MACHINE\SOFTWARE\Wow6432Node\Microsoft\VisualStudio\SxS\VS7
SET "VCINSTALLDIR=%ProgramFiles(x86)%\Microsoft Visual Studio 10.0\VC\"
SET "VSINSTALLDIR=%ProgramFiles(x86)%\Microsoft Visual Studio 10.0\"
:: if 64 bit
IF EXIST "%FrameworkDir64%" (
  SET "FrameworkDir=%FrameworkDir64%"
  SET "FrameworkVersion=%FrameworkVersion64%"
)

:Done_SetRegPath
:: --------------------------------------------------------------------------------------------
:: Set the CL.exe option to find the framework for linking
:: --------------------------------------------------------------------------------------------
IF "%CURRENT_CPU%"=="x86" (
  SET CL=/AI %FrameworkDir%%FrameworkVersion%
) ELSE (
  SET CL=/AI %FrameworkDir%\%FrameworkVersion%
)

:: --------------------------------------------------------------------------------------------
:: Save the default values for VCINSTALLDIR and VSINSTALLDIR just in case we get unexpected output from
:: the calls to REG in the next section
:: --------------------------------------------------------------------------------------------
SET "VCINSTALLDIR_Orig=%VCINSTALLDIR%"
SET "VSINSTALLDIR_Orig=%VSINSTALLDIR%"

:: --------------------------------------------------------------------------------------------
:: Read the value for VCINSTALLDIR and VSINSTALLDIR from the registry.
:: Note: The second call to REG will fail if VS is not installed.  These calls to REG are
:: checking to see if VS is installed in a custom location.  This behavior is expected in
:: this scenario.
:: --------------------------------------------------------------------------------------------
FOR /F "tokens=2* delims=	 " %%A IN ('REG QUERY "%RegKeyPath%" /v 10.0') DO SET VCINSTALLDIR=%%B
FOR /F "tokens=2* delims=	 " %%A IN ('REG QUERY "%VSRegKeyPath%" /v 10.0') DO SET VSINSTALLDIR=%%B

:: --------------------------------------------------------------------------------------------
:: Hide the error output from the call to 'REG' since VCINSTALLDIR and VSINSTALLDIR have already been set
:: to a default value.
:: --------------------------------------------------------------------------------------------
CLS

:: --------------------------------------------------------------------------------------------
:: Versions of Reg.exe on XP SP2 and SP3 have different output than other versions.  Detect
:: incorrect output and reset VCINSTALLDIR and VSINSTALLDIR as needed
:: --------------------------------------------------------------------------------------------
IF "%VCINSTALLDIR%"=="VERSION 3.0" SET "VCINSTALLDIR=%VCINSTALLDIR_Orig%"
IF "%VSINSTALLDIR%"=="VERSION 3.0" SET "VSINSTALLDIR=%VSINSTALLDIR_Orig%"

:: --------------------------------------------------------------------------------------------
:: Setup our VCTools environment path based on target CPU and processor architecture
:: When the native compilers are not installed for the specified CPU, attempt to locate the
:: cross-tools for non-native compilation.
:: --------------------------------------------------------------------------------------------
SET "VCTools=%VCINSTALLDIR%Bin"

IF "%CURRENT_CPU%" =="x64" (
  IF "%TARGET_CPU%" == "x64" (
    IF EXIST "%VCTools%\amd64\cl.exe" (
      SET "VCTools=%VCTools%\amd64;%VCTools%\VCPackages;"
      SET "SdkTools=%WindowsSdkDir%Bin\NETFX 4.0 Tools\x64;%WindowsSdkDir%Bin\x64;%WindowsSdkDir%Bin;"
      SET "FxTools=%FrameworkDir64%\%FrameworkVersion%;%FrameworkDir32%%FrameworkVersion%;%windir%\Microsoft.NET\Framework64\v3.5;%windir%\Microsoft.NET\Framework\v3.5;"
    ) ELSE (
      SET VCTools=
      ECHO The x64 compilers are not currently installed.
      ECHO Please go to Add/Remove Programs to update your installation.
      ECHO .
    )
  ) ELSE IF "%TARGET_CPU%" == "IA64" (
    IF EXIST "%VCTools%\x86_ia64\cl.exe" (
      SET "VCTools=%VCTools%\x86_ia64;%VCTools%\VCPackages;"
      SET "SdkTools=%WindowsSdkDir%Bin\NETFX 4.0 Tools\x64;%WindowsSdkDir%Bin\x64;%WindowsSdkDir%Bin;"
      SET "FxTools=%FrameworkDir64%\%FrameworkVersion%;%FrameworkDir32%%FrameworkVersion%;%windir%\Microsoft.NET\Framework64\v3.5;%windir%\Microsoft.NET\Framework\v3.5;"
    ) ELSE (
      SET VCTools=
      ECHO The IA64 compilers are not currently installed.
      ECHO Please go to Add/Remove Programs to update your installation.
      ECHO .
    )
  ) ELSE IF "%TARGET_CPU%" == "x86" (
    IF EXIST "%VCTools%\cl.exe" (
      SET "VCTools=%VCTools%;%VCTools%\VCPackages;"
      SET "SdkTools=%WindowsSdkDir%Bin\NETFX 4.0 Tools;%WindowsSdkDir%Bin;"
      SET "FxTools=%FrameworkDir32%%FrameworkVersion%;%windir%\Microsoft.NET\Framework\v3.5;"
    ) ELSE (
      SET VCTools=
      ECHO The x86 compilers are not currently installed.
      ECHO Please go to Add/Remove Programs to update your installation.
      ECHO .
    )
  )
) ELSE IF "%CURRENT_CPU%" =="IA64" (
  IF "%TARGET_CPU%" == "IA64" (
    IF EXIST "%VCTools%\IA64\cl.exe" (
      SET "VCTools=%VCTools%\IA64;%VCTools%;%VCTools%\VCPackages;"
      SET "SdkTools=%WindowsSdkDir%Bin\NETFX 4.0 Tools\IA64;%WindowsSdkDir%Bin\IA64;%WindowsSdkDir%Bin;"
      SET "FxTools=%FrameworkDir64%\%FrameworkVersion%;%FrameworkDir32%%FrameworkVersion%;%windir%\Microsoft.NET\Framework64\v3.5;%windir%\Microsoft.NET\Framework\v3.5;"
    ) ELSE (
      SET VCTools=
      ECHO The IA64 compilers are not currently installed.
      ECHO Please go to Add/Remove Programs to update your installation.
      ECHO .
    )
  ) ELSE IF "%TARGET_CPU%" == "x64" (
    IF EXIST "%VCTools%\x86_amd64\cl.exe" (
      SET "VCTools=%VCTools%\x86_amd64;%VCTools%\VCPackages;"
      SET "SdkTools=%WindowsSdkDir%Bin\NETFX 4.0 Tools\IA64;%WindowsSdkDir%Bin\IA64;%WindowsSdkDir%Bin;"
      SET "FxTools=%FrameworkDir64%\%FrameworkVersion%;%FrameworkDir32%%FrameworkVersion%;%windir%\Microsoft.NET\Framework64\v3.5;%windir%\Microsoft.NET\Framework\v3.5;"
    ) ELSE (
      SET VCTools=
      ECHO The VC compilers are not currently installed.
      ECHO Please go to Add/Remove Programs to update your installation.
      ECHO .
    )
  ) ELSE IF "%TARGET_CPU%" == "x86" (
    IF EXIST "%VCTools%\cl.exe" (
      SET "VCTools=%VCTools%;%VCTools%\VCPackages;"
      SET "SdkTools=%WindowsSdkDir%Bin\NETFX 4.0 Tools;%WindowsSdkDir%Bin;"
      SET "FxTools=%FrameworkDir32%%FrameworkVersion%;%windir%\Microsoft.NET\Framework\v3.5;"
    ) ELSE (
      SET VCTools=
      ECHO The x86 compilers are not currently installed.
      ECHO Please go to Add/Remove Programs to update your installation.
      ECHO .
    )
  )
) ELSE IF "%CURRENT_CPU%"=="x86" (
  IF "%TARGET_CPU%" == "x64" (
    IF EXIST "%VCTools%\x86_amd64\cl.exe" (
      SET "VCTools=%VCTools%\x86_amd64;%VCTools%\VCPackages;"
      SET "SdkTools=%WindowsSdkDir%Bin\NETFX 4.0 Tools;%WindowsSdkDir%Bin;"
      SET "FxTools=%FrameworkDir32%%FrameworkVersion%;%windir%\Microsoft.NET\Framework\v3.5;"
    ) ELSE (
      SET VCTools=
      ECHO The x64 cross compilers are not currently installed.
      ECHO Please go to Add/Remove Programs to update your installation.
      ECHO .
    )
  ) ELSE IF "%TARGET_CPU%" == "IA64" (
    IF EXIST "%VCTools%\x86_IA64\cl.exe" (
      SET "VCTools=%VCTools%\x86_IA64;%VCTools%;%VCTools%\VCPackages;"
      SET "SdkTools=%WindowsSdkDir%Bin\NETFX 4.0 Tools;%WindowsSdkDir%Bin;"
      SET "FxTools=%FrameworkDir32%%FrameworkVersion%;%windir%\Microsoft.NET\Framework\v3.5;"
    ) ELSE (
      SET VCTools=
      ECHO The IA64 compilers are not currently installed.
      ECHO Please go to Add/Remove Programs to update your installation.
      ECHO .
    )
  ) ELSE IF "%TARGET_CPU%" == "x86" (
    IF EXIST "%VCTools%\cl.exe" (
      SET "VCTools=%VCTools%;%VCTools%\VCPackages;"
      SET "SdkTools=%WindowsSdkDir%Bin\NETFX 4.0 Tools;%WindowsSdkDir%Bin;"
      SET "FxTools=%FrameworkDir32%%FrameworkVersion%;%windir%\Microsoft.NET\Framework\v3.5;"
    ) ELSE (
      SET VCTools=
      ECHO The x86 compilers are not currently installed. x86-x86
      ECHO Please go to Add/Remove Programs to update your installation.
      ECHO .
    )
  )
) ELSE IF EXIST "%VCTools%\cl.exe" (
  SET "VCTools=%VCTools%;%VCTools%\VCPackages;"
  SET "SdkTools=%WindowsSdkDir%Bin\NETFX 4.0 Tools;%WindowsSdkDir%Bin;"
  SET "FxTools=%FrameworkDir32%%FrameworkVersion%;%windir%\Microsoft.NET\Framework\v3.5;"
) ELSE (
  SET VCTools=
  ECHO The x86 compilers are not currently installed. default
  ECHO Please go to Add/Remove Programs to update your installation.
  ECHO .
)

:: --------------------------------------------------------------------------------------------
:: Set the values for VS directories.
:: --------------------------------------------------------------------------------------------
SET "DevEnvDir=%VSINSTALLDIR%Common7\IDE"
SET "VSTools=%VSINSTALLDIR%Common7\IDE;%VSINSTALLDIR%Common7\Tools;"
SET "VCLibraries=%VCINSTALLDIR%Lib"
SET "VCIncludes=%VCINSTALLDIR%INCLUDE"

:: --------------------------------------------------------------------------------------------
:: Display SDK environment information
:: --------------------------------------------------------------------------------------------
SET OSLibraries=%WindowsSdkDir%Lib
SET OSIncludes=%WindowsSdkDir%INCLUDE;%WindowsSdkDir%INCLUDE\gl

ECHO Setting SDK environment relative to %WindowsSdkDir%.
:: --------------------------------------------------------------------------------------------
:: Set whether this is a cross-compile enviornment
:: --------------------------------------------------------------------------------------------
IF "%CURRENT_CPU%"=="%TARGET_CPU%" (
  SET CommandPromptType=Native
) ELSE (
  SET CommandPromptType=Cross
)

SET Path=%FxTools%;%VSTools%;%VCTools%;%SdkTools%;%Path%

:: Clear the paths
SET LIB=
SET LIBPATH=
SET INCLUDE=

GOTO Set_%TARGET_CPU%

:Set_x86
IF NOT "x!VCINSTALLDIR!x"=="xx" (
  IF EXIST "!VCINSTALLDIR!ATLMFC" (
    SET "INCLUDE=!VCINSTALLDIR!ATLMFC\INCLUDE;%INCLUDE%"
    SET "LIB=!VCINSTALLDIR!ATLMFC\LIB;%LIB%"
  )
)

SET "LIB=%VCLibraries%;%OSLibraries%;%FxTools%;%LIB%"
SET "LIBPATH=%FxTools%;%VCLibraries%;%LIBPATH%"
SET "INCLUDE=%VCIncludes%;%OSIncludes%;%INCLUDE%"


GOTO Set_%TARGET_PLATFORM%

:Set_x64
IF NOT "x!VCINSTALLDIR!x"=="xx" (
  IF EXIST "!VCINSTALLDIR!ATLMFC\LIB\AMD64" (
    SET "INCLUDE=!VCINSTALLDIR!ATLMFC\INCLUDE;%INCLUDE%"
    SET "LIB=!VCINSTALLDIR!ATLMFC\LIB\AMD64;%LIB%"
  )
)

SET "LIB=%VCLibraries%\amd64;%OSLibraries%\X64;%LIB%"
SET "LIBPATH=%FxTools%;%VCLibraries%\amd64;%LIBPATH%"
SET "INCLUDE=%VCIncludes%;%OSIncludes%;%INCLUDE%"


GOTO Set_%TARGET_PLATFORM%

:Set_IA64
IF NOT "x!VCINSTALLDIR!x"=="xx" (
  IF EXIST "!VCINSTALLDIR!ATLMFC\LIB\IA64" (
    SET "INCLUDE=!VCINSTALLDIR!ATLMFC\INCLUDE;%INCLUDE%"
    SET "LIB=!VCINSTALLDIR!ATLMFC\LIB\IA64;%LIB%"
  )
)

SET "LIB=%VCLibraries%\IA64;%OSLibraries%\IA64;%LIB%"
SET "LIBPATH=%FxTools%;%VCLibraries%\ia64;%LIBPATH%"
SET "INCLUDE=%VCIncludes%;%OSIncludes%;%INCLUDE%"


GOTO Set_%TARGET_PLATFORM%

:: --------------------------------------------------------------------------------------------
:: Set Windows 7 specific variables
:: --------------------------------------------------------------------------------------------
:Set_WIN7
ECHO Targeting Windows 7 %TARGET_CPU% %DEBUGMSG%
ECHO.
SET APPVER=6.1
TITLE Microsoft Windows 7 %TARGET_CPU% %DEBUGMSG% Build Environment
GOTO End_Success

:: --------------------------------------------------------------------------------------------
:: Set Windows Server 2008 specific variables
:: --------------------------------------------------------------------------------------------
:Set_LHS
ECHO Targeting Windows Server 2008 %TARGET_CPU% %DEBUGMSG%
ECHO.
SET APPVER=6.0
TITLE Microsoft Windows Server 2008 %TARGET_CPU% %DEBUGMSG% Build Environment
GOTO End_Success

:: --------------------------------------------------------------------------------------------
:: Set Windows Vista  specific variables
:: --------------------------------------------------------------------------------------------
:Set_LH
ECHO Targeting Windows Vista %TARGET_CPU% %DEBUGMSG%
ECHO.
SET APPVER=6.0
TITLE Microsoft Windows Vista %TARGET_CPU% %DEBUGMSG% Build Environment
GOTO End_Success

:: --------------------------------------------------------------------------------------------
:: Set Windows Server 2003 specific variables
:: --------------------------------------------------------------------------------------------
:Set_SRV
ECHO Targeting Windows Server 2003 %TARGET_CPU% %DEBUGMSG%
ECHO.
SET APPVER=5.02
TITLE Microsoft Windows Server 2003 %TARGET_CPU% %DEBUGMSG% Build Environment
GOTO End_Success

:: --------------------------------------------------------------------------------------------
:: Set Windows XP specific variables
:: --------------------------------------------------------------------------------------------
:Set_XP
ECHO Targeting Windows XP %TARGET_CPU% %DEBUGMSG%
ECHO.
SET APPVER=5.01
TITLE Microsoft Windows XP %TARGET_CPU% %DEBUGMSG% Build Environment
GOTO End_Success

:: --------------------------------------------------------------------------------------------
:: Display command usage and goto cleanup code.
:: --------------------------------------------------------------------------------------------
:Error_Usage
ECHO Usage: "Setenv [/Debug | /Release][/x86 | /x64 | /ia64][/vista | /xp | /2003 | /2008 | /win7][-h | /?]"
ECHO.
ECHO                 /Debug   - Create a Debug configuration build environment
ECHO                 /Release - Create a Release configuration build environment
ECHO                 /x86     - Create 32-bit x86 applications
ECHO                 /x64     - Create 64-bit x64 applications
ECHO                 /ia64    - Create 64-bit ia64 applications
ECHO                 /vista   - Windows Vista applications
ECHO                 /xp      - Create Windows XP SP2 applications
ECHO                 /2003    - Create Windows Server 2003 applications
ECHO                 /2008    - Create Windows Server 2008 or Vista SP1 applications
ECHO                 /win7    - Create Windows 7 applications
ECHO.
ECHO Note:
ECHO * Platform(x86/x64/ia64) and PlatformToolSet(v90/v100/WindowsSDK7.1) set in project or solution will override the environment
ECHO * To upgrade VC6 or later projects to VC2010 format use the VCUpgrade.exe tool.
SET VCBUILD_DEFAULT_OPTIONS=
GOTO CleanUp

:: --------------------------------------------------------------------------------------------
:: End Successfully.  
:: If Windows 7 headers,libs and tools are not used, display a warning message.
:: If necessary, display warning about compiling on Windows 9x platforms.
:: --------------------------------------------------------------------------------------------
:End_Success

IF "x%OS%x" == "xWindows_NTx" SET "DEBUGMSG=" & GOTO CleanUp

ECHO *** WARNING ***
ECHO You are currently building on a Windows 9x based platform.  Most samples have 
ECHO NMAKE create a destination directory for created objects and executables.  
ECHO There is a known issue with the OS where NMAKE fails to create this destination
ECHO directory when the current directory is several directories deep.  To fix this 
ECHO problem, you must create the destination directory by hand from the command 
ECHO line before calling NMAKE. 
ECHO.

:DisplayWarningMsg_NoVersion
ECHO **********************************************************************************
ECHO WARNING: The VC++ Compiler Toolset is not using Windows SDK v7.1. 
ECHO **********************************************************************************
EXIT /B 0
:: -------------------------------------------------------------------
:: Persist Old Values and
:: -------------------------------------------------------------------
::SET TARGET_PLATFORM=
::SET Configuration=
::SET CURRENT_CPU=
::SET TARGET_CPU=
:: -------------------------------------------------------------------
:: Clean up
:: -------------------------------------------------------------------
:CleanUp
echo CURRENT_CPU = %CURRENT_CPU%
IF /I "%TARGET_CPU%"=="x86" GOTO x86
IF /I "%TARGET_CPU%"=="x64" GOTO x64

:x86
ECHO Building for x86
SET TCLDIR=tcltk
SET MACHINE=IX86
SET INSTALLDIR=./win32VC-tk8.6
GOTO run

:x64
ECHO Building for x64
SET TCLDIR=tcltk64
SET MACHINE=X64
SET INSTALLDIR=./win32VC-x86_64-tk8.6
GOTO run

:run
nmake -f makefile.vc MACHINE=%MACHINE% DEBUG=%DEBUG%
nmake -f makefile.vc install DEBUG=%DEBUG% INSTALLDIR=%INSTALLDIR%

SET OSLibraries=
SET OSIncludes=
SET VCINSTALLDIR_Orig=
SET VSINSTALLDIR_Orig=
SET VCTools=
SET VSTools=
SET FxTools=
SET SdkTools=
SET VCLibraries=
SET VCIncludes=
SET CurrentSdkVersion=
SET RegKeyPath=
SET VSRegKeyPath=
SET WinSDKRegKeyPath=
Set VCLibraries=
Set VCIncludes=
SET VCINSTALLDIR=
SET VSINSTALLDIR=
SET DevEnvDir=
SET FrameworkDir=
SET FrameworkDir64=
SET FrameworkVersion32=
SET FrameworkVersion64=
SET FrameworkDir32=
SET Framework35Version=
