import sys, os, re

for module in sys.argv[1: ]:
    py_dir = '/tmp/test_python_' + module
    os.system('rm -rf ' + py_dir)
    os.system(sys.executable + ' -m virtualenv ' + py_dir)
    os.system(py_dir + '/bin/pip install --no-cache-dir ' + module)
    os.system(py_dir + '/bin/python -m ' + module + '.test')
