from polynomial import Polynomial, Monomial
import solutionsToGroebnerBasis
import coordinates
import snappy

import re
import tempfile
import subprocess

class MagmaPrimaryIdeal(list):
    def __init__(self, polys, dimension = None, size = None):
        super(MagmaPrimaryIdeal,self).__init__(polys)
        self.dimension = dimension
        self.size = size

    def __repr__(self):
        return (
            "MagmaPrimaryIdeal([%s], dimension = %d, size = %d)" %
            (", ".join([repr(poly) for poly in self]),
            self.dimension,
            self.size))
        
def _eraseLineWraps(text):

    def processLineWithPotentialBackslash(line):
        strippedLine = line.strip()
        if strippedLine and strippedLine[-1] == '\\':
            return strippedLine[:-1]
        else:
            return strippedLine + '\n'

    return ''.join([processLineWithPotentialBackslash(line)
                    for line in text.split('\n')])

def parse_Magma(output):

    text = _eraseLineWraps(output)

    py_eval_match = re.search(
        r"PY=EVAL=SECTION=BEGINS=HERE(.*)PY=EVAL=SECTION=ENDS=HERE",
        text,
        re.DOTALL)

    assert py_eval_match, (
        "File not recognized as magma output (missing eval section)")
    py_eval = eval(py_eval_match.group(1))

    primary_decomposition_match = re.search(
        r"PRIMARY=DECOMPOSITION=BEGINS=HERE(.*)PRIMARY=DECOMPOSITION=ENDS=HERE",
        text,
        re.DOTALL)

    groebner_basis_match = re.search(
        r"GROEBNER=BASIS=BEGINS=HERE(.*)GROEBNER=BASIS=ENDS=HERE",
        text,
        re.DOTALL)

    if primary_decomposition_match:
    
        components_matches = re.findall(
            r"Ideal of Polynomial ring.*?"
            "Dimension (\d+).*?"
            "(Size of variety over algebraically closed field: (\d+).*?)?"
            "Groebner basis:\s*"
            "\[([^\]]*)\]",
            primary_decomposition_match.group(1),
            re.DOTALL)

        def parseInt(s):
            if s:
                return int(s)
            return None

        components = [
            MagmaPrimaryIdeal(
                polys = [ Polynomial.parseString(p)
                          for p in poly_strs.replace('\n',' ').split(',') ],
                dimension = int(dimension_str),
                size = parseInt(size_str))
            for dimension_str, variety_str, size_str, poly_strs
            in components_matches]
        
        return components, py_eval

    elif groebner_basis_match:

        polys_match = re.match(r"\s*\[([^\]]*)\]\s*", groebner_basis_match.group(1), re.DOTALL)
        assert polys_match
        polys_str = polys_match.group(1)

        polys = [ Polynomial.parseString(p)
                  for p in polys_str.replace('\n',' ').split(',') ]

        return [ MagmaPrimaryIdeal(polys) ], py_eval

    raise (
        "File not recognized as magma output "
        "(missing primary decomposition or groebner basis)")

def triangulation_from_Magma(output):
    """
    Reads the output from a magma computation and extracts the manifold for
    which this output constains solutions.
    """

    text = _eraseLineWraps(output)

    triangulation_match = re.search('==TRIANGULATION=BEGINS=='
                                    '(.*?)'
                                    '==TRIANGULATION=ENDS==',
                                    text,re.DOTALL)
    assert triangulation_match

    return snappy.Manifold(triangulation_match.group(1).strip())

def triangulation_from_Magma_file(filename):
    """
    Reads the output from a magma computation from the file with the given
    filename and extracts the manifold for which the file contains solutions.
    """

    return triangulation_from_Magma(open(filename).read())

def solutions_from_Magma_file(filename):

    """
    Reads the output from a magma computation from the file with the given
    filename and returns a list of solutions. Also see solutions_from_Magma.
    A non-zero dimensional component of the variety is reported as None.
    """

    return solutions_from_Magma(open(filename).read())

def solutions_from_Magma(output):
    """
    Assumes the given string is the output of a magma computation, parses
    it and returns a list of solutions.
    A non-zero dimensional component of the variety is reported as None.
    """

    components, extra_data = parse_Magma(output)

    def process_component(component):
        solutions = solutionsToGroebnerBasis.exact_solutions_with_one(
            component)
        
        if not component.dimension is None:
            if component.dimension > 0:
                assert solutions == [None]
        return solutions

    solutions = sum([process_component(component) for component in components],
                    [ ])

    def process_solution(solution):
        if not solution is None:
            return coordinates.PtolemyCoordinates(
                extra_data["variable_dict"](solution), is_numerical = False)
        return None

    solutions = [process_solution(solution)
                 for solution in solutions]

    return solutions

def run_Magma(content, filename_base, memory_limit, directory, verbose):

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
    out_file = resolved_dir + '/' + filename_base + '.magma_out'

    if verbose:
        print "Writing to file:", in_file

    open(in_file, 'w').write(content)

    if verbose:
        print "Magma's output in:", out_file

    cmd = 'ulimit -m %d; magma < "%s" > "%s"' % (
            int(memory_limit / 1024), in_file, out_file)

    if verbose:
        print "Command:", cmd
        print "Starting magma..."

    retcode = subprocess.call(cmd, shell = True)

    result = open(out_file, 'r').read()

    if verbose:
        print "magma finished."
        print "Parsing magma result..."

    return solutions_from_Magma(result)

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
