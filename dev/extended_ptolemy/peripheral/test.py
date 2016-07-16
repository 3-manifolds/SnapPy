import sys, getopt
import doctest

from . import dual_cellulation, link, peripheral, surface
modules = [dual_cellulation, link, peripheral, surface]

def verbose():
    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
        opts = [o[0] for o in optlist]
        verbose = '-v' in opts
    except getopt.GetoptError:
        verbose = False
    return verbose

def doctest_globals(module):
    if hasattr(module, 'doctest_globals'):
        return module.doctest_globals()
    else:
        return dict()

if __name__ == '__main__':
    failed, attempted = 0, 0
    for module in modules:
        print(module.__name__)
        result = doctest.testmod(module,
                                 extraglobs= doctest_globals(module),
                                 verbose=verbose())
        print(4*' ' + repr(result))
        failed += result.failed
        attempted += result.attempted
    print('\nAll doctests:\n    %s failures out of %s tests.' % (failed, attempted))

    
        
