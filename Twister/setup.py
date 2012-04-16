from distutils.core import setup, Extension

main_path = './lib/'
main_src = ['twister_coremodule.cpp']
kernel_path = './lib/kernel/'
kernal_src = ['twister.cpp', 'manifold.cpp', 'parsing.cpp', 'global.cpp']

core = Extension(
	name = 'twister.twister_core',
	sources = [main_path + file for file in main_src] + [kernel_path + file for file in kernal_src],
	include_dirs=[kernel_path],
	language='c++' )

setup(name='twister',
	version='2.3',
	description='Twister',
	author='Mark Bell',
	author_email='M.C.Bell@warwick.ac.uk',
	url='http://www.surfacebundles.wordpress.com/',
	packages=['twister'],
	package_dir={'twister':'lib'},
	ext_modules=[core] )
