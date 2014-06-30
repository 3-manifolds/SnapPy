import re

import snappy

"""
Basic functions to read a ptolemy solutions file.
"""

def erase_line_wraps(text):
    """
    If a line ends in \, the \ and \n is erased.
    """

    def process_line_with_potential_backslash(line):
        strippedLine = line.strip()
        if strippedLine and strippedLine[-1] == '\\':
            return strippedLine[:-1]
        else:
            return strippedLine + '\n'

    return ''.join([process_line_with_potential_backslash(line)
                    for line in text.split('\n')])

def find_section(text, name):

    """
    Finds all sections of the form
    
    NAME=BEGINS=HERE
    stuff
    stuff
    stuff
    NAME=ENDS=HERE
    
    or 

    ==NAME=BEGINS==
    stuff
    stuff
    stuff
    ==NAME=ENDS==

    in text where NAME has to be replaced by name.

    >>> t = (
    ... "abc abc\\n"
    ... "==FOO=BEGINS==\\n"
    ... "bar bar\\n"
    ... "==FOO=ENDS==\\n")
    >>> find_section(t, "FOO")
    ['bar bar']
    """

    old_style_regex = (name + "=BEGINS=HERE" +
                       "(.*?)" +
                       name + "=ENDS=HERE")
    new_style_regex = ("==" + name + "=BEGINS==" +
                       "(.*?)" +
                       "==" + name + "=ENDS==")
    regexs = [ old_style_regex, new_style_regex ]

    return [ s.strip()
             for regex in regexs
             for s in re.findall(regex, text, re.DOTALL) ]

def find_unique_section(text, name):

    """
    Like find_section but ensures that there is only a single
    such section. Exceptions in all other cases.

    >>> t = "abc abc\\n==FOO=BEGINS==\\nbar bar\\n==FOO=ENDS==\\n"
    >>> find_unique_section(t, "FOO")
    'bar bar'
    """

    sections = find_section(text, name)
    if len(sections) > 1:
        raise Exception("Section %s more than once in file" % name)
    if not sections:
        raise Exception("No section %s in file" % name)
    return sections[0]

def remove_outer_square_brackets(text):
    """
    Checks that test is of the form "[...]" and returns result between
    brackets.

    >>> remove_outer_square_brackets("[a*b*c]")
    'a*b*c'
    >>> remove_outer_square_brackets("[[a,b]")
    '[a,b'
    """

    if text[0] != '[' or text[-1] != ']':
        raise ValueError("Error while parsing: outer square brackets missing")

    return text[1:-1]

def parse_int_or_empty(s):

    """
    >>> parse_int_or_empty('3')
    3
    >>> parse_int_or_empty('') is None
    True
    """


    if s:
        return int(s)
    return None

def get_py_eval(text):

    """
    From a ptolemy solutions file, extract the PY=EVAL=SECTION
    """

    return eval(
        find_unique_section(erase_line_wraps(text), "PY=EVAL=SECTION"))

def get_manifold_thunk(text):

    """
    From a ptolemy solutions file, extract the manifold.
    Returned as thunk that evaluates to a snappy manifold.
    """

    def get_manifold(text = text):
        return snappy.Manifold(
            find_unique_section(erase_line_wraps(text), "TRIANGULATION"))

    return get_manifold

def get_manifold(text):

    """
    From a ptolemy solutions file, extract the manifold.
    Returned as snappy Manifold.
    """

    return get_manifold_thunk(text)()

def get_manifold_from_file(filename):

    """
    As get_manifold but takes filename.
    """

    return get_manifold(open(filename).read())
