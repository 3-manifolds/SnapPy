from distutils.core import setup, Extension, Command
import os

# So we can access all of the test suite just by doing "python setup.py test"
class TestCommand(Command):
	user_options = [ ]
	
	def initialize_options(self):
		pass
	
	def finalize_options(self):
		pass
	
	def run(self):
		''' Runs all of the test suite. '''
		from test.test import test_suite
		test_suite()

main_src = ['./lib/py_wrapper.cpp']
kernel_path = './lib/kernel/'
kernal_src = ['twister.cpp', 'manifold.cpp', 'parsing.cpp', 'global.cpp']

core = Extension(
	name = 'twister.twister_core',
	sources = main_src + [os.path.join(kernel_path, file) for file in kernal_src],
	include_dirs=[kernel_path],
	language='c++'
	)

setup(
	name='twister',
	version='2.4.0',
	description='Twister',
	author='Mark Bell',
	author_email='M.C.Bell@warwick.ac.uk',
	url='https://bitbucket.org/Mark_Bell/twister/',
	packages=['twister'],
	package_dir={'twister':'lib'},
	package_data={'twister': ['surfaces/*']},
	ext_modules=[core],
	cmdclass = {'test': TestCommand}
	)
