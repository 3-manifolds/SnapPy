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

import doctest, re, types

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

def cyopengl_replacement():
    """
    Have to run this late to avoid (circular?) import issues.
    """
    try:
        import snappy.CyOpenGL as CyOpenGL
        CYOPENGL = ''
    except ImportError:
        CYOPENGL = '#doctest: +SKIP'
    return CYOPENGL

if _within_sage:
    class DocTestParser(doctest.DocTestParser):
        def parse(self, string, name='<string>'):
            string = re.subn('#doctest: \+CYOPENGL', cyopengl_replacement(), string)[0]
            string = re.subn('([\n\A]\s*)sage:', '\g<1>>>>', string)[0]
            return doctest.DocTestParser.parse(self, string, name)

    globs = {'PSL':sage.all.PSL, 'BraidGroup':sage.all.BraidGroup}
else:
    class DocTestParser(doctest.DocTestParser):
        def parse(self, string, name='<string>'):
            string = re.subn('#doctest: \+CYOPENGL', cyopengl_replacement(), string)[0]
            return doctest.DocTestParser.parse(self, string, name)
        
    globs = dict()

def print_results(module, results):
    print(module.__name__ + ':')
    print('   %s failures out of %s tests.' %  (results.failed, results.attempted))
    
def doctest_modules(modules, verbose=False, print_info=True, extraglobs=dict()):
    finder = doctest.DocTestFinder(parser=DocTestParser())
    full_extraglobals = dict(globs.items() + extraglobs.items())
    failed, attempted = 0, 0
    for module in modules:
        if isinstance(module, types.ModuleType):
            runner = doctest.DocTestRunner(verbose=verbose)
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
    

        
    

