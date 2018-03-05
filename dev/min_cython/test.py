from distutils.core import run_setup

# Clean up
import os
os.system('rm -rf *.so build')

run_setup('setup.py', ['build_ext', '--inplace'])
import cymodule_one, cymodule_two

print('\n')

print(cymodule_one.timestwo(6))
A = cymodule_one.TwoByTwoMatrix([[1, 2],[3,4]])
B = cymodule_one.TwoByTwoMatrix([[5, 6],[7,8]])
print(A*B)
print(cymodule_two.add_matrices(A, B))
