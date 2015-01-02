import sys, os, re

if sys.argv[1].startswith('easy'):
    install_cmd = '/bin/easy_install '
elif sys.argv[1].startswith('pip'):
    install_cmd = '/bin/pip --no-cache-dir install '
else:
    print('Need to specify install via "pip" or "easy*"')
    sys.exit()


for module in sys.argv[2: ]:
    py_dir = '/tmp/test_python_' + module
    os.system('rm -rf ' + py_dir)
    os.system(sys.executable + ' -m virtualenv ' + py_dir)
    os.system(py_dir + install_cmd + module)
    os.system(py_dir + '/bin/python -m ' + module + '.test')
