def check(cmd, mf):
    m = mf.findNode('SnapPy')
    if m is None or m.filename is None:
        return None
    for pkg in [
        'os', 'sys', 'operator', 'types', 're', 'gzip',
        'struct', 'tempfile', 'tarfile', 'signal',  'ctypes', 'ctypes.util',
        'plink', 'polyviewer', 'horoviewer', 'manifolds',
        'numpy'
        ]:
        mf.import_hook(pkg, m, ['*'])
    return dict(
        packages = ['SnapPy']
    )
