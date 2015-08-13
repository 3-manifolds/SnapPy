import re

###########################################
### Methods to prepare and quote ASCII text

def _break_line_iterator(line, line_length):
    for i in range(0, len(line), line_length):
        yield line[i:i+line_length]

def break_long_lines(text, line_length = 76):
    """
    Break lines in ASCII text text longer than line_length by inserting 
    newline characters preceded by backslash (thus the resulting lines might
    be one character longer than ling_length).

    The reverse operation is join_long_lines, so
    join_long_lines(break_long_lines(text)) is always returning text as long as
    consisted of ASCII characters (even when text already had lines
    ending in a backslash as those backslashs are specially escaped)!

    This is consistent with the interpretation of the backslash newline
    sequence by many languages such as magma, python, C, C++ which treat these
    escaped newlines as non-existing. In particular, the result of
    break_long_lines and the text itself is the same to these languages (as
    long as the text didn't already contain lines ending in a backslash).

    >>> text = "This is a long line.\\nThis is an even longer line.\\n"
    >>> print(break_long_lines(text,8))
    This is \\
    a long l\\
    ine.
    This is \\
    an even \\
    longer l\\
    ine.
    <BLANKLINE>

    >>> join_long_lines(break_long_lines(text, 8)) == text
    True
    """

    def split_ending_backslash(line):
        if len(line) > 0 and line[-1] == '\\':
            return line[:-1], '\\\\\n'
        else:
            return line, ''

    def process_line(line):
        line_without, ending_backslash = split_ending_backslash(line)
        
        return ('\\\n'.join(_break_line_iterator(line_without, line_length)) +
                ending_backslash)

    return '\n'.join([process_line(line) for line in text.split('\n')])

def join_long_lines(text):
    """
    Deletes all backslash newline sequences. Inverse of break_long_lines.
    """

    return text.replace('\\\n','')

def join_long_lines_deleting_whitespace(text):
    """
    Similar to join_long_lines, but also deletes whitespace following a
    backslash newline sequence.
    Programs such as magma break long integers by introducing a backslash
    newline sequence and inserting extra whitespace for formatting.
    join_long_lines_deleting_whitespace can be used to undo this and parse
    the input normally.

    >>> join_long_lines_deleting_whitespace("Text:\\\\\\n   More")
    'Text:More'
    >>> join_long_lines_deleting_whitespace("    1234\\\\\\n    5678").strip()
    '12345678'
    """

    return re.sub(r'\\\n\s*','', text)


def quote_ascii_text(text):
    """
    Put the text in double quotes after escapes newlines, backslashes and
    double quotes. Giving the result of quote_ascii_text to eval should give
    the original string back if the string contained only ASCII characters.
    Similarly, giving the result of quote_ascii_text to magma's print,
    should give the original string back (magma's print might wrap long lines
    though).
    
    >>> text = 'Backslash:\\, Newline:\\n, Quote: "'
    >>> quote_ascii_text(text)
    '"Backslash:\\\\\\\\, Newline:\\\\n, Quote: \\\\""'
    >>> eval(quote_ascii_text(text)) == text
    True
    """

    def process_char(char):
        if char == '\n':
            return '\\n'
        if char == '\\':
            return '\\\\'
        if char == '"':
            return '\\"'
        return char

    return '"' + ''.join([process_char(c) for c in text]) + '"'

##########################################################################
### Iterators going through all tuples of integers with a certain property

def _lists_with_fixed_sum_iterator(N, l):
    if l == 1:
        yield [ N ]
    else:
        for i in range(N + 1):
            for j in _lists_with_fixed_sum_iterator(N-i, l-1):
                yield [i] + j

def tuples_with_fixed_sum_iterator(N, l, skipVertices = False):
    """
    Iterates through all l-tuples of non-negative integers summing up to N in
    lexicographic order. If skipVertices is True, N-tuples containing N, i.e.,
    of the form (0...0,N,0...0), are skipped.

    >>> list(tuples_with_fixed_sum_iterator(2, 3))
    [(0, 0, 2), (0, 1, 1), (0, 2, 0), (1, 0, 1), (1, 1, 0), (2, 0, 0)]
    >>> list(tuples_with_fixed_sum_iterator(2, 3, skipVertices = True))
    [(0, 1, 1), (1, 0, 1), (1, 1, 0)]

    """

    for i in _lists_with_fixed_sum_iterator(N, l):
        if not (skipVertices and (N in i)):
            yield tuple(i)

def triples_with_fixed_sum_iterator(N, skipVertices = False):
    """
    Similar to tuples_with_fixed_sum_iterator for triples.

    >>> list(triples_with_fixed_sum_iterator(2, skipVertices = True))
    [(0, 1, 1), (1, 0, 1), (1, 1, 0)]
    """

    return tuples_with_fixed_sum_iterator(N, 3, skipVertices = skipVertices)

def quadruples_with_fixed_sum_iterator(N, skipVertices = False):

    """
    Similar to tuples_with_fixed_sum_iterator for quadruples.

    >>> list(quadruples_with_fixed_sum_iterator(2, skipVertices = True))
    [(0, 0, 1, 1), (0, 1, 0, 1), (0, 1, 1, 0), (1, 0, 0, 1), (1, 0, 1, 0), (1, 1, 0, 0)]
    """

    return tuples_with_fixed_sum_iterator(N, 4, skipVertices = skipVertices)

############################################################################
### A subclass of Python list such that calling a method calls the method on
### every object in the list.

def _flatten(l, depth = 1):

    """
    Flatten a list or subclass of list. It will flatten the
    list until the given depth.
    >>> _flatten([0, [1, 2], [[3, 4], 5]], depth = 1)
    [0, 1, 2, [3, 4], 5]
    >>> _flatten([0, [1, 2], [[3, 4], 5]], depth = 2)
    [0, 1, 2, 3, 4, 5]
    """

    if depth == 0:
        return l

    result = []

    for e in l:
        if isinstance(e, list):
            result += _flatten(e, depth - 1)
        else:
            result.append(e)

    if isinstance(l, MethodMappingList):
        return type(l)(result, p = l)

    return result

class MethodMappingList(list):
    
    """
    Like a list but allows calling a method on it means that it is called
    for all its elements.

    >>> a = MethodMappingList([2+1j, 3+2j, 4+2j])
    >>> a.conjugate()
    [(2-1j), (3-2j), (4-2j)]

    This can be nested:
    
    >>> b = MethodMappingList([1+1j, a])
    >>> b.conjugate()
    [(1-1j), [(2-1j), (3-2j), (4-2j)]]

    Also supports flattening:

    >>> b.flatten()
    [(1+1j), (2+1j), (3+2j), (4+2j)]
    
    """

    def __init__(self, l = [], p = None):
        super(MethodMappingList, self).__init__(l)

    def __call__(self, *args, **kwargs):
        
        return type(self)([elt(*args, **kwargs) for elt in self],
                          p = self)

    def __getattr__(self, attr):

        return type(self)([getattr(e, attr) for e in self],
                          p = self)

    def flatten(self, depth = 1):
        return _flatten(self, depth = depth)


def _test():
    """
    >>> text = "Line\\\\nLong Line\\\\\\nVery long line and...\\nHere \\n\\n \\n\\n \\n \\\\n"
    >>> for i in range(1, 80):
    ...     if not join_long_lines(break_long_lines(text, i)) == text:
    ...         raise Exception("Test failure")
    """

    pass
    
