# -*- mode: python -*-
from PyInstaller.utils.hooks import collect_submodules, collect_data_files
import sys
cpu_width = '64bit' if sys.maxsize > 2**32 else '32bit'

block_cipher = None

options = [('v', None, 'OPTION')]

imports = collect_submodules('snappy')
imports += collect_submodules('cypari')
imports += collect_submodules('jedi')

datafiles = collect_data_files('jedi')
datafiles += collect_data_files('parso')
datafiles += collect_data_files('snappy_manifolds')
datafiles += collect_data_files('snappy'
datafiles += collect_data_files('plink')
datafiles += collect_data_files('spherogram')


# SnapPyHP.pyd and twister_core.pyd are compiled with the MS Visual
# C++ compiler from Visual Studio 2015, which dynamically links them
# against the C++ runtime library msvcp140.dll.  As of PyInstaller
# 3.3, vcruntime140.dll (and msvcp100.dll) are listed in
# PyInstaller.depends.dylib._includes but msvcp140.dll is not.  To
# work around this we add msvcp140.dll as a binary, specifying
# that it should be at the top level of the bundle, adjacent to the
# two .pyd files which depend on it.

if cpu_width == '32bit':
    binaries = [('C:\\Windows\\SysWOW64\\msvcp140.dll', '.')]
else:
    binaries = [('C:\\Windows\\System32\\msvcp140.dll', '.')]

a = Analysis(['SnapPy.py'],
             binaries=binaries,
             hiddenimports=imports + ['linecache'],
             datas=datafiles,
             hookspath=[],
             runtime_hooks=[],
             excludes=['gi', 'pytz', 'td', 'sphinx', 'alabaster', 'babel',
                       'idlelib', 'bsddb'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)

# As of PyInstaller 3.3, the system dlls api-ms-win-core*.dll and
# api-ms-win-crt*.dll are explicity included in the app bundle (see
# pyinstaller.depend.dylib._includes).  But these DLLs are not
# compatible across different versions of Windows, so an app built on
# Windows 10 will crash on Windows 7.  Moreover, they do not need to
# be included in the app bundle because they are always available from
# the OS.  Until this is fixed, the following hack is a workaround.

a.binaries = [b for b in a.binaries if b[1].find('api-ms-win') < 0]

pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='SnapPy',
          debug=True,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='SnapPy_debug')
