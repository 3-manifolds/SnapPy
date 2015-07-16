"""
Speeding up compilation of SnapPy for development purposes.  Just do

python dev/ccache.py command [options]

as you would with "setup.py".  For now, assumes Unix-like system with
gcc (not clang) as the compiler and relies on having "ccache" installed.

Unfortunately, no dependencies are deduced or checked here.  Do
"ccache -C" to clear the cache.  In particular, if you edit something
in "kernel" the corresponding "quad_double" code will not be
recompiled unless you manually touch said file.

---------

I tried to use the parallel compilation code given at:

http://stackoverflow.com/questions/11013851/speeding-up-build-process-with-distutils

without any success.
"""
import os, sys
import multiprocessing
import distutils.unixccompiler
import distutils.spawn
import distutils.errors 

# Make use ccache 
os.environ['CC'] = 'ccache gcc'
os.environ['CXX'] = 'ccache gcc++'

# monkey-patch for parallel compilation

def compile_one(command):
    try:
        distutils.spawn.spawn(command)
    except distutils.errors.DistutilsExecError, msg:
        return msg

def compile_parallel(self, sources, output_dir=None, macros=None, include_dirs=None,
                     debug=0, extra_preargs=None, extra_postargs=None, depends=None):
    # Next two lines from distutils.ccompiler.CCompiler.compile:
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(
            output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)

    # Next two lines from distutils.unixccompiler.UnixCCompiler._compile
    compiler_so = self.compiler_so
    if sys.platform == 'darwin':
        compiler_so = distutils.unixccompiler._osx_support.compiler_fixup(
            compiler_so, cc_args + extra_postargs)
    
    compile_commands = []
    for obj in objects:
        if obj in build:
            src, ext = build[obj]
            compile_commands.append(compiler_so + cc_args + [src, '-o', obj] + extra_postargs)

    #num_proc = 2
    #pool = multiprocessing.Pool(num_proc)
    #map = pool.map 
    for msg in map(compile_one, compile_commands):
        if not msg is None:
            raise distutils.errors.CompileError(msg)

    # Return *all* object filenames, not just the ones we just built.
    return objects


#distutils.unixccompiler.UnixCCompiler.compile = compile_parallel

if __name__ == '__main__':
    # Run the usual setup
    sys.path = [os.path.abspath(os.curdir)] + sys.path
    import setup
