docs = """
To install a package from PyPI via pip and run it's test suite:

    python test_pypi.py -p module

To use easy install instead, give the "-e" flag; to prefer the
versions of packages on TestPyPI add "-t".   Typically examples:

   python test_pypi.py -p -t snappy
   python test_pypi.py -e -t plink spherogram snappy
"""
   

import sys, os, re, subprocess, argparse
import virtualenv, setuptools   # just to make sure these are installed

parser = argparse.ArgumentParser(description='Check packages on (Test)PyPI via virtualenvs.',
                                 epilog=docs, formatter_class=argparse.RawDescriptionHelpFormatter)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-p', '--pip', help='Use pip to install the packages',
                    action='store_true')
group.add_argument('-e', '--easy_install', help='Use easy_install to aquire the packages',
                    action='store_true')
parser.add_argument('-s', '--source', help='Make pip not use wheels',
                    action='store_true')
parser.add_argument('-t', '--testing', help='Use testingpypi not real pypi', action='store_true')
parser.add_argument('modules', nargs='+')


class Sandbox:
    def __init__(self, name):
        if sys.platform.startswith('win'):
            tmp_dir, bin_dir, exe = 'tmp', 'Scripts', '.exe'
        else:
            tmp_dir, bin_dir, exe = '/tmp', 'bin', ''

        py_dir = os.path.join(tmp_dir, 'test_python_' + name)
        if os.path.exists(py_dir):
            print('Deleting existing virtualenv')
            os.system('rm -rf ' + py_dir)

        print('Creating virtualenv in ' + py_dir)
        os.system(sys.executable + ' -m virtualenv ' + py_dir)
        self.bin_dir, self.py_dir, self.exe = bin_dir, py_dir, exe

    def execute(self, command):
        command[0] = os.path.join(self.py_dir, self.bin_dir, command[0] + self.exe)
        subprocess.call(command)

if __name__ == '__main__':
    args = parser.parse_args()
    testpypi = 'https://testpypi.python.org/simple'
    if args.pip:
        install_cmd = ['pip', 'install', '--no-cache-dir', '--pre']
        if args.testing:
            install_cmd += ['--extra-index-url', testpypi]
        if args.source:
            install_cmd += ['--no-use-wheel']
    elif args.easy_install:
        install_cmd = ['easy_install']


    sandbox = Sandbox(args.modules[-1])
    for module in args.modules:
        if args.easy_install and args.testing:
            sandbox.execute(install_cmd + ['-f', testpypi + '/' + module, module])
        else:
            sandbox.execute(install_cmd + [module])

    for module in args.modules:
        sandbox.execute(['python', '-m', module + '.test'])


