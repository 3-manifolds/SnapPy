docs = """
To install a package from PyPI via pip and run it's test suite:

    python test_pypi.py module

To prefer the versions of packages on TestPyPI add "-t"; to ignore
wheels and build from source and "-s".  Typical examples:

   python test_pypi.py -t FXrays
   python test_pypi.py -s snappy
"""


import sys, os, re, shutil, subprocess, argparse
import setuptools   # just to make sure these are installed
if not sys.platform.startswith('linux'):
    import venv

parser = argparse.ArgumentParser(description='Check packages on (Test)PyPI via venvs.',
                                 epilog=docs, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-s', '--source', help='Make pip not use wheels',
                    action='store_true')
parser.add_argument('-t', '--testing', help='Use testingpypi not real pypi', action='store_true')

if sys.platform == 'darwin':
    parser.add_argument('-u', '--universal2', help='Favor universal2 wheels if available', action='store_true')

parser.add_argument('-w', '--wheelhouse', help='Local directory of wheels', default='')

parser.add_argument('modules', nargs='+')

# cf. https://bugs.python.org/issue22490
environ = os.environ.copy()
environ.pop('__PYVENV_LAUNCHER__', None)


class Sandbox:
    def __init__(self, name):
        if sys.platform.startswith('win'):
            tmp_dir, bin_dir, exe = 'tmp', 'Scripts', '.exe'
        else:
            tmp_dir, bin_dir, exe = '/tmp', 'bin', ''

        sys.version_info
        py_dir = os.path.join(tmp_dir,
                              'test_python_%d.%d.%d_' % sys.version_info[:3] + name)
        if os.path.exists(py_dir):
            print('Deleting existing venv')
            shutil.rmtree(py_dir)

        print('Creating venv in ' + py_dir)
        if sys.platform.startswith('linux'):
            subprocess.call(['venv', '-p', sys.executable, py_dir], env=environ)
        else:
            subprocess.call([sys.executable, '-m', 'venv', py_dir], env=environ)

        self.bin_dir, self.py_dir, self.exe = bin_dir, py_dir, exe
        self.site_packages = os.path.join(py_dir, 'lib',
                                          'python%d.%d' % sys.version_info[:2],
                                          'site-packages')

    def execute(self, command):
        command[0] = os.path.join(self.py_dir, self.bin_dir, command[0] + self.exe)
        subprocess.call(command, env=environ)

if __name__ == '__main__':
    args = parser.parse_args()
    testpypi = 'https://test.pypi.org/simple'
    install_cmd = ['pip', 'install', '--no-cache-dir', '--pre']
    sandbox = Sandbox(args.modules[-1])
    if args.wheelhouse:
        install_cmd += ['--find-links=' + args.wheelhouse]
    if args.testing:
        install_cmd += ['--extra-index-url', testpypi]
    if args.source:
        install_cmd += ['--no-binary=' + ','.join(args.modules)]
    else:
        if sys.platform=='darwin' and args.universal2:
            install_cmd += ['--platform=macosx_10_9_universal2',
                            '--only-binary=:all:',
                            '--target=' + sandbox.site_packages]

    for module in args.modules:
        sandbox.execute(install_cmd + [module])

    for module in args.modules:
        sandbox.execute(['python', '-m', module + '.test'])
