"""
Helper code for dealing with additional functionality when Sage is
present.

Any method which works only in Sage should be decorated with
"@sage_method" and any doctests (in Sage methods or not) which should
be run only in Sage should be styled with input prompt "sage:" rather
than the usual ">>>".
"""

try:
    import sage.all
    _within_sage = True
except:
    _within_sage = False
    import decorator

import sys, doctest, re, types

from .numeric_output_checker import NumericOutputChecker

class SageNotAvailable(Exception):
    pass

if _within_sage:
    def sage_method(function):
        function._sage_method = True
        return function


    try: # Sage >= 9.3, see https://trac.sagemath.org/ticket/24483
        from sage.rings.complex_mpfr import (ComplexField,
                                             ComplexField_class,
                                             create_ComplexNumber)
    except ModuleNotFoundError:
        from sage.rings.complex_field import ComplexField, ComplexField_class
        from sage.rings.complex_number import create_ComplexNumber

else:
    def _sage_method(function, *args, **kw):
        raise SageNotAvailable('Sorry, this feature requires using SnapPy inside Sage.')

    def sage_method(function):
        return decorator.decorator(_sage_method, function)


# Not currently used, but could be exploited by an interpreter to hide
# sage_methods when in plain Python.

def sage_methods(obj):
    ans = []
    for attr in dir(obj):
        try:
            methods = getattr(obj, attr)
            if methods._sage_method is True:
                ans.append(methods)
        except AttributeError:
            pass
    return ans

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
    _use_cyopengl_initialized = False
    _use_cyopengl = False
    use_modernopengl = True
    use_sage = False

    def parse(self, string, name='<string>'):
        string = re.subn(
            r'#doctest: \+CYOPENGL',
            '' if DocTestParser._use_cyopengl else '#doctest: +SKIP',
            string)[0]

        string = re.subn(
            r'#doctest: \+CYMODERNOPENGL',
            (''
             if (DocTestParser._use_cyopengl and
                 DocTestParser.use_modernopengl)
             else '#doctest: +SKIP'),
            string)[0]
        
        if DocTestParser.use_sage:
            string = re.subn(r'(\n\s*)sage:|(\A\s*)sage:',
                             r'\g<1>>>>',
                             string)[0]
        return doctest.DocTestParser.parse(self, string, name)

DocTestParser.use_sage = _within_sage

if _within_sage:
    globs = {'PSL':sage.all.PSL, 'BraidGroup':sage.all.BraidGroup}
else:
    globs = { }

def print_results(module, results):
    root = tk_root()
    # Platform specific hacks to make running the tests work.
    if root and (sys.version_info.major < 3 or not root_is_fake()):
        if sys.platform in ('linux', 'linux2'):
            root.deiconify()
            root.update_idletasks()
        else:
            root.update()
    print(module.__name__ + ':')
    print('   %s failures out of %s tests.' %  (results.failed, results.attempted))

def doctest_modules(modules, verbose=False, print_info=True, extraglobs=dict()):
    finder = doctest.DocTestFinder(parser=DocTestParser())
    #full_extraglobals = dict(globs.items() + extraglobs.items())
    full_extraglobals = globs.copy()
    full_extraglobals.update(extraglobs)
    failed, attempted = 0, 0
    for module in modules:
        if isinstance(module, types.ModuleType):
            runner = doctest.DocTestRunner(checker = NumericOutputChecker(), verbose=verbose)
            for test in finder.find(module, extraglobs=full_extraglobals):
                runner.run(test)
            result = runner.summarize()
        else:
            result = module(verbose=verbose)
        failed += result.failed
        attempted += result.attempted
        if print_info:
            print_results(module, result)

    if print_info:
        print('\nAll doctests:\n   %s failures out of %s tests.' % (failed, attempted))
    return doctest.TestResults(failed, attempted)
