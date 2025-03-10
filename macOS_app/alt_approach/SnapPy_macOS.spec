# -*- mode: python -*-
from PyInstaller.utils.hooks import collect_submodules, collect_data_files
import sys

block_cipher = None

options = [('v', None, 'OPTION')]

imports = collect_submodules('snappy')
imports += collect_submodules('snappy_15_knots')
imports += collect_submodules('cypari')
imports += collect_submodules('jedi')
imports += collect_submodules('pyx')
imports += collect_submodules('low_index')
imports += collect_submodules('tkinter_gl')


datafiles = collect_data_files('jedi')
datafiles += collect_data_files('pyx')
datafiles += collect_data_files('parso')
datafiles += collect_data_files('snappy_manifolds')
datafiles += collect_data_files('snappy_15_knots')
datafiles += collect_data_files('snappy')
datafiles += collect_data_files('plink')
datafiles += collect_data_files('spherogram')
datafiles += collect_data_files('tkinter_gl')


# SnapPyHP.pyd and twister_core.pyd are compiled with the MS Visual
# C++ compiler from Visual Studio 2015, which dynamically links them
# against the C++ runtime library msvcp140.dll.  As of PyInstaller
# 3.3, vcruntime140.dll (and msvcp100.dll) are listed in
# PyInstaller.depends.dylib._includes but msvcp140.dll is not.  To
# work around this we add msvcp140.dll as a binary, specifying
# that it should be at the top level of the bundle, adjacent to the
# two .pyd files which depend on it.

a = Analysis(['SnapPy.py'],
             hiddenimports=imports + ['linecache', 'pkg_resources.py2_warn', 'pkg_resources.extern'],
             datas=datafiles,
             hookspath=[],
             runtime_hooks=[],
             excludes=['gi', 'pytz', 'td', 'sphinx', 'alabaster', 'babel',
                       'idlelib', 'bsddb'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)

# As of PyInstaller 3.3, the system dlls api-ms-win-core*.dll and
# api-ms-win-crt*.dll are explicitly included in the app bundle (see
# pyinstaller.depend.dylib._includes).  But these DLLs are not
# compatible across different versions of Windows, so an app built on
# Windows 10 will crash on Windows 7.  Moreover, they do not need to
# be included in the app bundle because they are always available from
# the OS.  Until this is fixed, the following hack is a workaround.

# a.binaries = [b for b in a.binaries if b[1].find('api-ms-win') < 0]

pyz = PYZ(a.pure,
          a.zipped_data,
          cipher=block_cipher)

exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='SnapPy',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=False,
          disable_windowed_traceback=False,
          target_arch='universal2',
          codesign_identity=None,
          entitlements_file=None,
          icon='../icons/SnapPy.icns')

coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas, 
               strip=False,
               upx=True,
               upx_exclude=[],
               name='SnapPy')

app = BUNDLE(coll,
             name='SnapPy.app',
             icon='../icons/SnapPy.icns',
             bundle_identifier=None)
