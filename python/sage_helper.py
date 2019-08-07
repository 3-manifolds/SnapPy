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

from .numericOutputChecker import NumericOutputChecker

class SageNotAvailable(Exception):
    pass

if _within_sage:
    def sage_method(function):
        function._sage_method = True
        return function
else:
    def _sage_method(function, *args, **kw):
        raise SageNotAvailable('Sorry, this feature requires using SnapPy inside Sage.')

    def sage_method(function):
        return decorator.decorator(_sage_method, function)


# Not currently used, but could be exploited by an interpeter to hide
# sage_methods when in plain Python.

def sage_methods(obj):
    ans = []
    for attr in dir(obj):
        try:
            methods = getattr(obj, attr)
            if methods._sage_method == True:
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
    # if we are running the tests from the snappy app the default root will
    # already exist -- it will be the tkterminal window.  Otherwise, we open
    # a root window here to serve as the master of all of the GUI windows
    # which get created during testing.
    if _gui_status['tk'] and not Tk_._default_root:
        try:
            root = Tk_.Tk()
            root.wm_geometry('-100+100')
            root.wait_visibility()
            Tk_.Label(root, text='Close me when done.').pack(padx=20, pady=20)
            root.update_idletasks()
            _gui_status['fake_root'] = True
            if sys.version_info.major < 3:
                Tk_._default_root = root
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

if _within_sage:
    class DocTestParser(doctest.DocTestParser):
        def parse(self, string, name='<string>'):
            if not hasattr(self, 'cyopengl_replacement'):
                self.cyopengl_replacement = '' if cyopengl_works() else '#doctest: +SKIP'
            string = re.subn('#doctest: \+CYOPENGL', self.cyopengl_replacement, string)[0]
            string = re.subn('(\n\s*)sage:|(\A\s*)sage:', '\g<1>>>>', string)[0]
            return doctest.DocTestParser.parse(self, string, name)

    globs = {'PSL':sage.all.PSL, 'BraidGroup':sage.all.BraidGroup}
else:
    class DocTestParser(doctest.DocTestParser):
        def parse(self, string, name='<string>'):
            if not hasattr(self, 'cyopengl_replacement'):
                self.cyopengl_replacement = '' if cyopengl_works() else '#doctest: +SKIP'
            string = re.subn('#doctest: \+CYOPENGL', self.cyopengl_replacement, string)[0]
            return doctest.DocTestParser.parse(self, string, name)

    globs = dict()

def print_results(module, results):
    root = tk_root()
    # Hack to mitigate hangs when running the tests from the app in linux.
    if root and sys.version_info.major < 3 or not root_is_fake():
        root.update_idletasks()
        root.deiconify()
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
