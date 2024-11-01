from .numeric_output_checker import NumericOutputChecker
from .sage_helper import _within_sage

import doctest
import getopt
import re
import sys
import types

# Used for doctesting
_gui_status = {}

try:
    from snappy.gui import Tk_
    _gui_status['tk'] = True
except ImportError:
    _gui_status['tk'] = False
if _gui_status['tk']:
    try:
        import snappy.CyOpenGL
        _gui_status['cyopengl'] = True
    except:
        _gui_status['cyopengl'] = False
else:
    _gui_status['cyopengl'] = False
_gui_status['fake_root'] = False

def cyopengl_works():
    if not _gui_status['cyopengl']:
        return False
    # if we are running the tests from the snappy app the default root will
    # already exist -- it will be the tkterminal window.  Otherwise, we open
    # a root window here to serve as the master of all of the GUI windows
    # which get created during testing.
    if _gui_status['tk'] and not Tk_._default_root:
        try:
            root = Tk_.Tk()
            if sys.platform not in ('linux', 'linux2'):
                root.withdraw()
        except:
            # tkinter loads OK but is not able to get a display.
            _gui_status['tk'] = _gui_status['cyopengl'] = False
    return _gui_status['cyopengl']


def tk_root():
    if _gui_status['tk']:
        return Tk_._default_root
    else:
        return None

def root_is_fake():
    return _gui_status['fake_root']

class DocTestParser(doctest.DocTestParser):
    use_cyopengl = False
    use_cymodernopengl = True
    use_sage = _within_sage

    def parse(self, string, name='<string>'):
        string = re.subn(
            r'#doctest: \+CYOPENGL',
            '' if DocTestParser.use_cyopengl else '#doctest: +SKIP',
            string)[0]

        string = re.subn(
            r'#doctest: \+CYMODERNOPENGL',
            (''
             if (DocTestParser.use_cyopengl and
                 DocTestParser.use_cymodernopengl)
             else '#doctest: +SKIP'),
            string)[0]

        if DocTestParser.use_sage:
            string = re.subn(r'(\n\s*)sage:|(\A\s*)sage:',
                             r'\g<1>>>>',
                             string)[0]
        return doctest.DocTestParser.parse(self, string, name)

if _within_sage:
    try:
        from sage.all import PSL, BraidGroup
    except ImportError:
        import sage.groups.perm_gps.permgroup_element
        from sage.groups.perm_gps.permgroup_named import PSL
        from sage.groups.braid import BraidGroup
    globs = {'PSL': PSL, 'BraidGroup': BraidGroup}
else:
    globs = {}

def print_results(module, results):
    root = tk_root()
    # Platform specific hacks to make running the tests work.
    if root and not root_is_fake():
        if sys.platform in ('linux', 'linux2'):
            root.deiconify()
            root.update_idletasks()
        else:
            root.update()
    print(module.__name__ + ':')
    print('   %s failures out of %s tests.' % (results.failed,
                                               results.attempted))

def doctest_modules(modules, verbose=False, print_info=False, extraglobs={}):
    finder = doctest.DocTestFinder(parser=DocTestParser())
    full_extraglobals = globs.copy()
    full_extraglobals.update(extraglobs)
    failed, attempted = 0, 0
    for module in modules:
        if isinstance(module, types.ModuleType):
            runner = doctest.DocTestRunner(checker=NumericOutputChecker(), verbose=verbose)
            for test in finder.find(module, extraglobs=full_extraglobals):
                runner.run(test)
            result = runner.summarize()
        else:
            result = module(verbose=verbose, print_info=False)
        failed += result.failed
        attempted += result.attempted
        if print_info:
            print_results(module, result)

    if print_info:
        print('\nAll doctests:\n   %s failures out of %s tests.' % (failed, attempted))
    return doctest.TestResults(failed, attempted)

def run_doctests_as_main(run_doctests):
    optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
    verbose = len(optlist) > 0
    results = run_doctests(verbose=verbose, print_info=True)
    sys.exit(results.failed)

