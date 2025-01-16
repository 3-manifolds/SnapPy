#! /usr/bin/env python3

"""
Adds doc directory to exisiting wheel
"""

import subprocess
import sys
import os
import re
import tempfile
import glob
import shutil

doc_zipfile = os.path.abspath(sys.argv[1])
wheel_names = glob.glob(sys.argv[2] + '/snappy-*.whl')
if len(wheel_names) == 0:
    raise ValueError('No snappy wheels in wheeldir!')

def normalize(path):
    return set(re.split(r'\-|\.', os.path.basename(path)))

for wheel_path in wheel_names:
    python = sys.executable
    tmp_dir = tempfile.mkdtemp()
    subprocess.check_call([python, '-m', 'wheel', 'unpack', '--dest', tmp_dir, wheel_path])
    pkg, version = os.path.basename(wheel_path).split('-')[:2]
    wheel_dir = os.path.join(tmp_dir, pkg + '-' + version)
    target = os.path.join(wheel_dir, pkg)
    target_doc_dir = os.path.join(target, 'doc')
    if os.path.exists(target_doc_dir):
        print('Deleting existing docs...')
        shutil.rmtree(target_doc_dir)
    print('Unpacking docs..')
    subprocess.check_call(['unzip', '-q', doc_zipfile, '-d', target])
    subprocess.check_call([python, '-m', 'wheel', 'pack', '--dest', tmp_dir, wheel_dir])
    new_whl = glob.glob(os.path.join(tmp_dir,  pkg + '-' + version + '*.whl'))[0]
    if normalize(new_whl) != normalize(wheel_path):
        raise ValueError(f'Wheel tag changed, see {tmp_dir}')
    subprocess.check_call(['cp', new_whl, wheel_path])
    print(f'Added docs to {os.path.basename(wheel_path)} in-place')
    



        
