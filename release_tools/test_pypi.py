import sys, os, re, subprocess

if sys.argv[1].startswith('easy'):
    cmd = ['easy_install']
elif sys.argv[1].startswith('pip'):
    cmd = ['pip', 'install', '--no-cache-dir']
else:
    print('Need to specify install via "pip" or "easy*"')
    sys.exit()


for module in sys.argv[2: ]:
    if sys.platform.startswith('win'):
        tmp_dir, bin_dir, exe = 'tmp', 'Scripts', '.exe'
    else:
        tmp_dir, bin_dir = '/tmp', 'bin', ''

    py_dir = os.path.join(tmp_dir, 'test_python_' + module)
    os.system('rm -rf ' + py_dir)
    os.system(sys.executable + ' -m virtualenv ' + py_dir)
    cmd[0] = os.path.join(py_dir, bin_dir, cmd[0] + exe)
    cmd.append(module)
    subprocess.call(cmd)
    test = [os.path.join(py_dir, bin_dir, 'python' + exe), 
            '-m', module + '.test']
    subprocess.call(test)
