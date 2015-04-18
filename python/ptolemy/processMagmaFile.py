from __future__ import print_function
from .polynomial import Polynomial
from .ptolemyVarietyPrimeIdealGroebnerBasis import PtolemyVarietyPrimeIdealGroebnerBasis
from . import processFileBase
from . import utilities
import snappy

import re
import tempfile
import subprocess
        
###############################################################################
# functions

def decomposition_from_magma(text):
    
    py_eval = processFileBase.get_py_eval(text)
    manifold_thunk = processFileBase.get_manifold_thunk(text)

    untyped_decomposition = processFileBase.find_section(
        text, "IDEAL=DECOMPOSITION")
    primary_decomposition = processFileBase.find_section(
        text, "PRIMARY=DECOMPOSITION")
    radical_decomposition = processFileBase.find_section(
        text,  "RADICAL=DECOMPOSITION")

    if untyped_decomposition:
        decomposition = untyped_decomposition[0]
    elif primary_decomposition:
        decomposition = primary_decomposition[0]
    elif radical_decomposition:
        decomposition = radical_decomposition[0]
    else:
        raise ValueError(
            "File not recognized as magma output "
            "(missing primary decomposition or radical decomposition)")
    
    # Remove "\" at the wrapped lines
    decomposition = utilities.join_long_lines_deleting_whitespace(
        decomposition)
    
    # Remove outer square brackets
    decomposition = processFileBase.remove_outer_square_brackets(decomposition)

    decomposition_comps = [ c.strip() for c in decomposition.split(']') ]
    decomposition_components = [ c + ']' for c in decomposition_comps if c ]

    free_variables_section = processFileBase.find_section(
        text, "FREE=VARIABLES=IN=COMPONENTS")
    if free_variables_section:
        free_variables = eval(free_variables_section[0])
    else:
        free_variables = len(decomposition_components) * [ None ]

    witnesses_section = processFileBase.find_section(
        text, "WITNESSES=FOR=COMPONENTS")
    if witnesses_section:
        witnesses_sections = processFileBase.find_section(
            witnesses_section[0], "WITNESSES")
    else:
        witnesses_sections = len(decomposition_components) * [ "" ]

    def process_match(i, comp, free_vars, witnesses_txt):

        if i != 0:
            if not comp[0] == ',':
                raise ValueError("Parsing decomposition, expected "
                                 "separating comma.")
            comp = comp[1:].strip()

        witnesses_txts = processFileBase.find_section(
            witnesses_txt, "WITNESS")
        
        witnesses = [
            _parse_ideal_groebner_basis(
                utilities.join_long_lines_deleting_whitespace(t).strip(),
                py_eval, manifold_thunk, free_vars, [])
            for t in witnesses_txts ]

        return _parse_ideal_groebner_basis(comp, py_eval, manifold_thunk,
                                           free_vars, witnesses)
        
    return utilities.MethodMappingList(
        [ process_match(i, comp, free_vars, witnesses)
          for i, (comp, free_vars, witnesses)
          in enumerate(zip(decomposition_components,
                           free_variables, witnesses_sections)) ])

def _parse_ideal_groebner_basis(text, py_eval, manifold_thunk,
                                free_vars, witnesses):
    match = re.match(
        r"Ideal of Polynomial ring of rank.*?\n"
        "\s*?(Order:\s*?(.*?)|(.*?)\s*?Order)\n"
        "\s*?Variables:(.*?\n)+"
        ".*?Dimension (\d+).*?\s*([^,]*[Pp]rime)?.*?\n"
        "(\s*?Size of variety over algebraically closed field: (\d+).*?\n)?"
        "\s*Groebner basis:\n"
        "\s*?\[([^\[\]]*)\]$",
        text)
        
    if not match:
        raise ValueError("Parsing error in component of "
                         "decomposition: %s" % text)

    (
        tot_order_str, post_order_str, pre_order_str,
        var_str,
        dimension_str, prime_str,
        variety_str, size_str,
        
        poly_strs                                      ) = match.groups()
        
    dimension = int(dimension_str)
        
    if dimension == 0:
        polys = [ Polynomial.parse_string(p)
                  for p in poly_strs.replace('\n',' ').split(',') ]
    else:
        polys = []
            
    order_str = post_order_str if post_order_str else pre_order_str
    if not order_str:
        raise ValueError("Could not parse order in decomposition")

    if order_str.strip().lower() == 'lexicographical':
        term_order = 'lex'
    else:
        term_order = 'other'

    is_prime = (prime_str is None) or (prime_str.lower() == 'prime')

    return  PtolemyVarietyPrimeIdealGroebnerBasis(
        polys = polys,
        term_order = term_order,
        size = processFileBase.parse_int_or_empty(size_str),
        dimension = dimension,
        is_prime = is_prime,
        free_variables = free_vars,
        py_eval = py_eval,
        manifold_thunk = manifold_thunk,
        witnesses = witnesses)


def triangulation_from_magma(text):
    """
    Reads the output from a magma computation and extracts the manifold for
    which this output constains solutions.
    """

    return processFileBase.get_manifold(text)

def triangulation_from_magma_file(filename):
    """
    Reads the output from a magma computation from the file with the given
    filename and extracts the manifold for which the file contains solutions.
    """

    return processFileBase.get_manifold_from_file(filename)

def contains_magma_output(text):
    return ("IDEAL=DECOMPOSITION=BEGINS" in text or
            "PRIMARY=DECOMPOSITION=BEGINS" in text or
            "RADICAL=DECOMPOSITION=BEGINS" in text)

def solutions_from_magma_file(filename, numerical = False):

    """
    Obsolete, use processFileDispatch.parse_solutions_from_file instead.

    Reads the output from a magma computation from the file with the given
    filename and returns a list of solutions. Also see solutions_from_magma.
    A non-zero dimensional component of the variety is reported as
    NonZeroDimensionalComponent.
    """

    return solutions_from_magma(open(filename).read(), numerical)

def solutions_from_magma(output, numerical = False):
    """
    Obsolete, use processFileDispatch.parse_solutions instead.

    Assumes the given string is the output of a magma computation, parses
    it and returns a list of solutions.
    A non-zero dimensional component of the variety is reported as
    NonZeroDimensionalComponent.
    """

    return decomposition_from_magma(output).solutions(
        numerical = numerical)

def run_magma(content,
              filename_base, memory_limit, directory, verbose):

    """
    call magma on the given content and 
    """

    if directory:
        resolved_dir = directory
        if not resolved_dir[-1] == '/':
            resolved_dir = resolved_dir + '/'
    else:
        resolved_dir = tempfile.mkdtemp() + '/'

    in_file  = resolved_dir + filename_base + '.magma'
    out_file = resolved_dir + filename_base + '.magma_out'

    if verbose:
        print("Writing to file:", in_file)

    open(in_file, 'w').write(content)

    if verbose:
        print("Magma's output in:", out_file)

    cmd = 'ulimit -m %d; echo | magma "%s" > "%s"' % (
            int(memory_limit / 1024), in_file, out_file)

    if verbose:
        print("Command:", cmd)
        print("Starting magma...")

    retcode = subprocess.call(cmd, shell = True)

    result = open(out_file, 'r').read()

    if verbose:
        print("magma finished.")
        print("Parsing magma result...")

    return decomposition_from_magma(result)

###############################################################################
# magma test

_magma_output_for_4_1__sl3 = """
==TRIANGULATION=BEGINS==
% Triangulation
4_1
geometric_solution  2.02988321
oriented_manifold
CS_known -0.0000000000000001

1 0
    torus   0.000000000000   0.000000000000

2
   1    1    1    1 
 0132 1302 1023 2031
   0    0    0    0 
  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  0  1  0 -1 -1  0  2 -1  0 -1  0  1  0  1 -1  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  0.500000000000   0.866025403784

   0    0    0    0 
 0132 1302 1023 2031
   0    0    0    0 
  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  0  1 -2  1  1  0  0 -1  0  1  0 -1  0 -1  1  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  0.500000000000   0.866025403784


==TRIANGULATION=ENDS==
PY=EVAL=SECTION=BEGINS=HERE
{'variable_dict' : 
     (lambda d, negation = (lambda x:-x): {
          'c_1020_0' : d['c_0012_1'],
          'c_1020_1' : d['1'],
          'c_0201_0' : d['c_0201_0'],
          'c_0201_1' : d['c_0201_0'],
          'c_2100_0' : d['1'],
          'c_2100_1' : d['c_0012_1'],
          'c_2010_0' : d['1'],
          'c_2010_1' : d['c_0012_1'],
          'c_0102_0' : d['c_0102_0'],
          'c_0102_1' : d['c_0102_0'],
          'c_1101_0' : d['c_1101_0'],
          'c_1101_1' : negation(d['c_1101_0']),
          'c_1200_0' : d['c_0012_1'],
          'c_1200_1' : d['1'],
          'c_1110_0' : negation(d['c_1011_1']),
          'c_1110_1' : negation(d['c_1011_0']),
          'c_0120_0' : d['c_0102_0'],
          'c_0120_1' : d['c_0102_0'],
          'c_2001_0' : d['c_0201_0'],
          'c_2001_1' : d['c_0201_0'],
          'c_0012_0' : d['1'],
          'c_0012_1' : d['c_0012_1'],
          'c_0111_0' : d['1'],
          'c_0111_1' : negation(d['1']),
          'c_0210_0' : d['c_0201_0'],
          'c_0210_1' : d['c_0201_0'],
          'c_1002_0' : d['c_0102_0'],
          'c_1002_1' : d['c_0102_0'],
          'c_1011_0' : d['c_1011_0'],
          'c_1011_1' : d['c_1011_1'],
          'c_0021_0' : d['c_0012_1'],
          'c_0021_1' : d['1']})}
PY=EVAL=SECTION=ENDS=HERE
PRIMARY=DECOMPOSITION=BEGINS=HERE
[
    Ideal of Polynomial ring of rank 7 over Rational Field
    Lexicographical Order
    Variables: t, c_0012_1, c_0102_0, c_0201_0, c_1011_0, c_1011_1, c_1101_0
    Dimension 0, Radical, Prime
    Size of variety over algebraically closed field: 2
    Groebner basis:
    [
        t - 3/8*c_1011_1 - 1/2,
        c_0012_1 + 1/2*c_1011_1 + 3/2,
        c_0102_0 + 1/2*c_1011_1 + 1/2,
        c_0201_0 + 1/2*c_1011_1 + 1/2,
        c_1011_0 - c_1011_1 - 3,
        c_1011_1^2 + 3*c_1011_1 + 4,
        c_1101_0 - 1
    ],
    Ideal of Polynomial ring of rank 7 over Rational Field
    Lexicographical Order
    Variables: t, c_0012_1, c_0102_0, c_0201_0, c_1011_0, c_1011_1, c_1101_0
    Dimension 0, Radical, Prime
    Size of variety over algebraically closed field: 2
    Groebner basis:
    [
        t - 1/2*c_1101_0 - 15/8,
        c_0012_1 - 1,
        c_0102_0 + 4/3*c_1101_0 - 2/3,
        c_0201_0 - 4/3*c_1101_0 - 1/3,
        c_1011_0 - 1/3*c_1101_0 - 1/3,
        c_1011_1 + 1/3*c_1101_0 + 1/3,
        c_1101_0^2 - 1/4*c_1101_0 + 1
    ],
    Ideal of Polynomial ring of rank 7 over Rational Field
    Lexicographical Order
    Variables: t, c_0012_1, c_0102_0, c_0201_0, c_1011_0, c_1011_1, c_1101_0
    Dimension 0, Radical, Prime
    Size of variety over algebraically closed field: 2
    Groebner basis:
    [
        t - c_1011_1 - 1,
        c_0012_1 - 1,
        c_0102_0 - c_1011_1,
        c_0201_0 - c_1011_1,
        c_1011_0 + c_1011_1,
        c_1011_1^2 + c_1011_1 + 1,
        c_1101_0 - 1
    ]
]
PRIMARY=DECOMPOSITION=ENDS=HERE
CPUTIME : 0.020

Total time: 0.419 seconds, Total memory usage: 5.62MB
"""
