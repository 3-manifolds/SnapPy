from __future__ import print_function
from .polynomial import Polynomial, Monomial
from .solutionsToGroebnerBasis import exact_solutions_with_one
from .numericalSolutionsToGroebnerBasis import numerical_solutions_with_one
from .component import NonZeroDimensionalComponent, ZeroDimensionalComponent
from . import coordinates
import snappy

try:
    from sage.libs.pari import gen 
    from sage.libs.pari.gen import pari
    from sage.rings.complex_field import ComplexField
    _within_sage = True
except ImportError:
    from cypari import gen
    from cypari.gen import pari
    _within_sage = False

import re
import tempfile
import subprocess

MAGMA_PRINT_ADDITIONAL_DATA = """

print "==TRIANGULATION" cat "=BEGINS==";
print "$QUOTED_TRIANGULATION";
print "==TRIANGULATION" cat "=ENDS==";
print "PY=EVAL=SECTION" cat "=BEGINS=HERE";
print "$PY_EVAL_SECTION";
print "PY=EVAL=SECTION=ENDS=HERE";

"""

MAGMA_PRIMARY_DECOMPOSITION_TEMPLATE = """
// Setting up the Polynomial ring and ideal

R<$VARIABLES_WITH_NON_ZERO_CONDITION> :=\
 PolynomialRing(RationalField(), $VARIABLE_WITH_NON_ZERO_CONDITION_NUMBER);
MyIdeal := ideal<R |
          $EQUATIONS_WITH_NON_ZERO_CONDITION>;

""" + MAGMA_PRINT_ADDITIONAL_DATA + """

// Value indicating failure
P := -1;

// Computing the primary decomposition
P,Q := PrimaryDecomposition(MyIdeal);

if Type(P) eq RngIntElt then
    // Some error occured
    print "PRIMARY=DECOMPOSITION" cat "=FAILED";
    exit;
else
    // Success
    print "PRIMARY=DECOMPOSITION" cat "=BEGINS=HERE";
    P;
    print "PRIMARY=DECOMPOSITION" cat "=ENDS=HERE";
end if;


"""

MAGMA_GROEBNER_BASIS_TEMPLATE = """
R<$VARIABLES> := PolynomialRing(RationalField(), $VARIABLE_NUMBER);
MyIdeal := ideal<R |
          $EQUATIONS>;

""" + MAGMA_PRINT_ADDITIONAL_DATA + """


// Value indicating failure
G := -1;

// Computing the Groebner basis
G := GroebnerBasis(MyIdeal);

if Type(G) eq RngIntElt then
    // Some error occured
    print "GROEBNER=BASIS" cat "=FAILED";
else
    // Success
    print "GROEBNER=BASIS" cat "=BEGINS=HERE";
    G;
    print "GROEBNER=BASIS" cat "=ENDS=HERE";
end if;


"""

MAGMA_VARIETY_TEMPLATE = """

// Setting up the Polynomial ring and ideal

R<$VARIABLES_WITH_NON_ZERO_CONDITION> :=\
 PolynomialRing(RationalField(), $VARIABLE_WITH_NON_ZERO_CONDITION_NUMBER);
MyIdeal := ideal<R |
          $EQUATIONS_WITH_NON_ZERO_CONDITION>;

""" + MAGMA_PRINT_ADDITIONAL_DATA + """

// Value indicating failure
G := -1;
P := -1;

totTime := Cputime();

if false then
    // Computing the Groebner basis
    groebTime := Cputime();
    G := GroebnerBasis(MyIdeal);
    print "GROEBNER_BASIS_TIME: ", Cputime(groebTime);
    if Type(G) eq RngIntElt then
        // Some Error 
        print "GROEBNER=BASIS" cat "=FAILED";
        print "TOTAL_TIME: ", Cputime(totTime);
        exit;
    else
        // Success
        print "GROEBNER=BASIS" cat "=BEGINS=HERE";
        G;
        print "GROEBNER=BASIS" cat "=ENDS=HERE";
    end if;
end if;

if false then
    // Computing the Radical decomposition
    radTime := Cputime();
    P := RadicalDecomposition(MyIdeal);
    print "RADICAL_DECOMPOSITION_TIME: ", Cputime(radTime);

    if Type(P) eq RngIntElt then
        // Some Error 
        print "RADICAL=DECOMPOSITION" cat "=FAILED";
        print "TOTAL_TIME: ", Cputime(totTime);
        exit;
    else
        // Success
        print "RADICAL=DECOMPOSITION" cat "=BEGINS=HERE";
        P;
        print "RADICAL=DECOMPOSITION" cat "=ENDS=HERE";
    end if;
else
    // Computing the primary decomposition
    primTime := Cputime();
    P, Q := PrimaryDecomposition(MyIdeal);
    print "PRIMARY_DECOMPOSITION_TIME: ", Cputime(primTime);

    if Type(P) eq RngIntElt then
        // Some Error 
        print "PRIMARY=DECOMPOSITION" cat "=FAILED";
        print "TOTAL_TIME: ", Cputime(totTime);
        exit;
    else
        // Success
        print "PRIMARY=DECOMPOSITION" cat "=BEGINS=HERE";
        P;
        print "PRIMARY=DECOMPOSITION" cat "=ENDS=HERE";
    end if;
end if;

// Decimal digits precision
precision := 100;

C<I> := ComplexField(precision);

print "PRECISION" cat "=BEGINS=HERE";
precision;
print "PRECISION" cat "=ENDS=HERE";
print "VARIABLE=ORDER" cat "=BEGINS=HERE";
print "$VARIABLES_WITH_NON_ZERO_CONDITION";
print "VARIABLE=ORDER" cat "=ENDS=HERE";

print "VARIETY" cat "=BEGINS=HERE";

varTime := Cputime();

VARIETY_FAILED := 0;

for Comp in P do
    D := Dimension(Comp);

    print "VARIETY=COMPONENT" cat "=BEGINS=HERE";

    if D eq 0 then
        V := -1;
        V := Variety(Comp, C);

        if Type(V) eq RngIntElt then
            VARIETY_FAILED := 1;
        else
            V;
        end if;

    else
        print "NonZeroDimensionalComponent(dimension = " * 
              IntegerToString(D) * ")";
    end if;

    print "VARIETY=COMPONENT" cat "=ENDS=HERE";

end for;

if VARIETY_FAILED eq 0 then
    print "VARIETY" cat "=ENDS=HERE";
end if;

print "PRIMARY_DECOMPOSITION_TIME: ", Cputime(varTime);
print "TOTAL_TIME: ", Cputime(totTime);

"""


class MagmaIdeal(list):
    def __init__(self, polys, dimension = None, size = None, primary = False):
        super(MagmaIdeal,self).__init__(polys)
        self.dimension = dimension
        self.size = size
        self.primary = primary

    def __repr__(self):
        return (
            "MagmaIdeal([%s], dimension = %d, size = %s, primary = %s)" %
            (", ".join([repr(poly) for poly in self]),
            self.dimension,
            self.size,
            self.primary))
        
def _eraseLineWraps(text):

    def process_line_with_potential_backslash(line):
        strippedLine = line.strip()
        if strippedLine and strippedLine[-1] == '\\':
            return strippedLine[:-1]
        else:
            return strippedLine + '\n'

    return ''.join([process_line_with_potential_backslash(line)
                    for line in text.split('\n')])

def _find_magma_section(text, section_name):

    """
    Finds a section of the form
    
    SECTION=NAME=BEGINS=HERE
    stuff
    stuff
    stuff
    SECTION=NAME=ENDS=HERE
    
    in text
    """

    return [
        s.strip()
        for s in re.findall(
            section_name + "=BEGINS=HERE(.*?)"  +section_name + "=ENDS=HERE",
            text, re.DOTALL)]

def _parse_magma_ideal(text):

    primary_decomposition = _find_magma_section(text,
                                                "PRIMARY=DECOMPOSITION")
    
    radical_decomposition = _find_magma_section(text,
                                                "RADICAL=DECOMPOSITION")

    groebner_basis = _find_magma_section(text,
                                         "GROEBNER=BASIS")

    if primary_decomposition or radical_decomposition:
    
        if primary_decomposition:
            decomposition = primary_decomposition[0]
        else:
            decomposition = radical_decomposition[0]

        def find_first_square_bracket_group(text):
            assert text[0] == '['
            nested = 0
            for i in range(0, len(text)):
                if text[i] == '[':
                    nested += +1
                if text[i] == ']':
                    nested += -1
                if nested == 0:
                    return text[:i+1]

            raise ValueError("Parsing Error")

        primary_decomposition_string = find_first_square_bracket_group(
            decomposition)

        components_matches = re.findall(
            r"Ideal of Polynomial ring.*?"
            "Dimension (\d+).*?"
            "(Size of variety over algebraically closed field: (\d+).*?)?"
            "Groebner basis:\s*"
            "\[([^\]]*?)\]",
            primary_decomposition_string,
            re.DOTALL)

        def parse_int(s):
            if s:
                return int(s)
            return None

        components = [
            MagmaIdeal(
                polys = [ Polynomial.parse_string(p)
                          for p in poly_strs.replace('\n',' ').split(',') ],
                dimension = int(dimension_str),
                size = parse_int(size_str),
                primary = True)
            for dimension_str, variety_str, size_str, poly_strs
            in components_matches]
        
        return components

    elif groebner_basis:

        polys_match = re.match(r"\s*\[([^\]]*)\]\s*",
                               groebner_basis[0], re.DOTALL)
        assert polys_match
        polys_str = polys_match.group(1)

        polys = [ Polynomial.parse_string(p)
                  for p in polys_str.replace('\n',' ').split(',') ]

        return [ MagmaIdeal(polys) ]

    raise ValueError(
        "File not recognized as magma output "
        "(missing primary decomposition or groebner basis)")

def triangulation_from_magma(output):
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

def triangulation_from_magma_file(filename):
    """
    Reads the output from a magma computation from the file with the given
    filename and extracts the manifold for which the file contains solutions.
    """

    return triangulation_from_magma(open(filename).read())

def solutions_from_magma_file(filename, numerical = False):

    """
    Reads the output from a magma computation from the file with the given
    filename and returns a list of solutions. Also see solutions_from_magma.
    A non-zero dimensional component of the variety is reported as
    NonZeroDimensionalComponent.
    """

    return solutions_from_magma(open(filename).read(), numerical)

def solutions_from_magma(output, numerical = False):
    """
    Assumes the given string is the output of a magma computation, parses
    it and returns a list of solutions.
    A non-zero dimensional component of the variety is reported as
    NonZeroDimensionalComponent.
    """

    text = _eraseLineWraps(output)

    py_eval_sections = _find_magma_section(text, "PY=EVAL=SECTION")
    assert len(py_eval_sections) == 1, (
        "File not recognized as magma output (missing eval section)")
    py_eval = eval(py_eval_sections[0])

    manifoldThunk = lambda : triangulation_from_magma(output)

    if numerical:

        variety_section = _find_magma_section(text, "VARIETY")
        variable_order = _find_magma_section(text, "VARIABLE=ORDER")

        if variety_section and variable_order:

            precision = int(_find_magma_section(text, "PRECISION")[0])

            return _get_numerical_solutions(variety_section[0],
                                            variable_order[0],
                                            precision,
                                            py_eval, manifoldThunk)

    components = _parse_magma_ideal(text)

    def process_component(component):
        if component.primary:
            if not component.dimension is None:
                if component.dimension > 0:
                    return [ NonZeroDimensionalComponent(
                            dimension = component.dimension) ]

        if numerical:
            raw_solutions = numerical_solutions_with_one(component)

            if component.primary:
                solutions = [ ZeroDimensionalComponent(raw_solutions) ]
            else:
                solutions = [ raw_solutions ]

        else:
            solutions = exact_solutions_with_one(component)
        
        if not component.dimension is None:
            if component.dimension > 0:
                assert len(solutions) == 1
                assert isinstance(solutions[0], NonZeroDimensionalComponent)
        return solutions

    solutions = sum([process_component(component) for component in components],
                    [ ])

    def process_solution(solution):
        if isinstance(solution, NonZeroDimensionalComponent):
            return solution

        def toPtolemyCoordinates(sol):
            return coordinates.PtolemyCoordinates(
                sol,
                is_numerical = numerical,
                py_eval_section = py_eval,
                manifoldThunk = manifoldThunk)

        if isinstance(solution, list):
            l = [ toPtolemyCoordinates(sol) for sol in solution ]
            if isinstance(solution, ZeroDimensionalComponent):
                return ZeroDimensionalComponent(l)
            return l
        return toPtolemyCoordinates(solution)
        
    solutions = [process_solution(solution)
                 for solution in solutions]
    
    return solutions

def run_magma(content,
              filename_base, numerical, memory_limit, directory, verbose):

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

    return solutions_from_magma(result, numerical = numerical)

def _get_numerical_solutions(variety_section,
                             variable_order,
                             precision,
                             py_eval,
                             manifoldThunk):

    if precision < pari.get_real_precision() - 2:
        raise ValueError("More precision requested than stored in magma "
                         "file. Reduce with pari.set_real_precision(...).")

    variables = [ var.strip() for var in variable_order.split(',') ]

    def process_component(component):
        
        nonZeroMatch = re.match(
            r'NonZeroDimensionalComponent\(\s*dimension\s*=\s*(\d+)\s*\)',
            component)

        if nonZeroMatch:
            return NonZeroDimensionalComponent(dimension = int(
                    nonZeroMatch.group(1)))

        assert component[0] == '['
        assert component[-1] == ']'
        
        points = re.findall('<(.*?)>', component, re.DOTALL)

        def process_point(point):
            
            numbers = [ pari(number) for number in point.split(',') ]

            d = dict(zip(variables, numbers))
            d['1'] = pari('1')

            return coordinates.PtolemyCoordinates(
                d, is_numerical = True, py_eval_section = py_eval,
                manifoldThunk = manifoldThunk)

        return [ process_point(point) for point in points ]

    return [
        process_component(component)
        for component in _find_magma_section(variety_section,
                                             "VARIETY=COMPONENT")]

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
