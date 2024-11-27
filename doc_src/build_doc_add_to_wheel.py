"""
Builds the docs with the current Python and adds them to the given wheel.
"""

import subprocess
import sys
import os
import tempfile
import glob
import shutil

wheel_path = os.path.abspath(sys.argv[1])

if not wheel_path.endswith('.whl'):
    raise ValueError('Usage: python build_doc_add_to_wheel.py snappy.whl')

python = sys.executable
tmp_dir = tempfile.mkdtemp()
subprocess.check_call([python, '-m', 'sphinx', '-b', 'html', '-E',
                       '-d', tmp_dir + '/doctrees', '.', tmp_dir + '/doc'])
subprocess.check_call([python, '-m', 'wheel', 'unpack', '--dest', tmp_dir, wheel_path])
pkg, version = os.path.basename(wheel_path).split('-')[:2]
wheel_dir = os.path.join(tmp_dir, pkg + '-' + version)
target = os.path.join(wheel_dir, pkg)
target_doc_dir = os.path.join(target, 'doc')
if os.path.exists(target_doc_dir):
    print('Deleting existing docs...')
    shutil.rmtree(target_doc_dir)
print('Moving docs.')
shutil.move(tmp_dir + '/doc', target_doc_dir)
subprocess.check_call([python, '-m', 'wheel', 'pack', '--dest', tmp_dir, wheel_dir])
new_whl = glob.glob(os.path.join(tmp_dir,  pkg + '-' + version + '*.whl'))[0]
if os.path.basename(new_whl) != os.path.basename(wheel_path):
    raise ValueError(f'Wheel tag changed, see {tmp_dir}')
subprocess.check_call(['cp', new_whl, wheel_path])
print(f'Added docs to {os.path.basename(wheel_path)} in-place')



        
