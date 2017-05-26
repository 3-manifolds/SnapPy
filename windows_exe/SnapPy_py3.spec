# -*- mode: python -*-
from PyInstaller.utils.hooks import collect_submodules, collect_data_files

block_cipher = None

options = [('v', None, 'OPTION')]

imports = collect_submodules('snappy')
imports += collect_submodules('cypari')
imports += collect_submodules('jedi')

datafiles = collect_data_files('jedi')

a = Analysis(['SnapPy.py'],
             binaries=None,
             hiddenimports=imports + ['linecache'],
             datas=datafiles,
             hookspath=[],
             runtime_hooks=[],
             excludes=['gi', 'pytz', 'td', 'sphinx', 'alabaster', 'babel',
                       'idlelib', 'bsddb'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)

# On line 99 of PyInstaller/building/build_main.py it says:
# "On Windows files from C:\Windows are excluded by default."
# But this is false when running in Python 3.6.  Until this is
# fixed, the following hack makes it true.  If these DLLs are
# included, apps built on Windows 10 will crash on Windows 7.

a.binaries = [b for b in a.binaries if b[1].find('api-ms-win') < 0]

pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='SnapPy',
          debug=False,
          strip=False,
          upx=True,
          console=False )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='SnapPy')
