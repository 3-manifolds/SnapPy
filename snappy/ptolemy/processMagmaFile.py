from __future__ import print_function
from .polynomial import Polynomial
from .ptolemyVarietyPrimeIdealGroebnerBasis import PtolemyVarietyPrimeIdealGroebnerBasis
from .component import MethodForwardingList
from . import processFileBase
import snappy

import re
import tempfile
import subprocess
        
###############################################################################
# functions

def decomposition_from_magma(text):
    
    py_eval = processFileBase.get_py_eval(text)
    manifold_thunk = processFileBase.get_manifold_thunk(text)

    text = processFileBase.erase_line_wraps(text)

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

    # Remove outer square brackets
    decomposition = processFileBase.remove_outer_square_brackets(decomposition)

    decomposition_comps = [ c.strip() for c in decomposition.split(']') ]
    decomposition_components = [ c for c in decomposition_comps if c ]

    free_variables_section = processFileBase.find_section(
        text, "FREE=VARIABLES=IN=COMPONENTS")
    if free_variables_section:
        free_variables = eval(free_variables_section[0])
    else:
        free_variables = len(decomposition_components) * [ None ]

    def process_match(i, comp, free_vars):

        if i != 0:
            if not comp[0] == ',':
                raise ValueError("Parsing decomposition, expected "
                                 "separating comma.")
            comp = comp[1:].strip()

        match = re.match(
            r"Ideal of Polynomial ring of rank.*?\n"
            "\s*?(Order:\s*?(.*?)|(.*?)\s*?Order)\n"
            "\s*?Variables:(.*?\n)+"
            ".*?Dimension (\d+).*?,\s*([^,]*[Pp]rime).*?\n"
            "(\s*?Size of variety over algebraically closed field: (\d+).*?\n)?"
            "\s*Groebner basis:\n"
            "\s*?\[([^\[\]]*)$",
            comp)
        
        if not match:
            raise ValueError("Parsing error in component of "
                             "decomposition: %s" % comp)

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

        is_prime = (prime_str.lower() == 'prime')

        return  PtolemyVarietyPrimeIdealGroebnerBasis(
            polys = polys,
            term_order = term_order,
            size = processFileBase.parse_int_or_empty(size_str),
            dimension = dimension,
            is_prime = is_prime,
            free_variables = free_vars,
            py_eval = py_eval,
            manifold_thunk = manifold_thunk)
        
    return MethodForwardingList(
        [ process_match(i, comp, free_vars)
          for i, (comp, free_vars)
          in enumerate(zip(decomposition_components, free_variables)) ])


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
# MAGMA template section

########################################
# MAGMA pieces

# MAGMA: print out the data needed besides the ideal to recover the Ptolemy
# coordinates and the triangulation
_MAGMA_PRINT_ADDITIONAL_DATA = """

print "==TRIANGULATION" cat "=BEGINS==";
print "$QUOTED_TRIANGULATION";
print "==TRIANGULATION" cat "=ENDS==";
print "PY=EVAL=SECTION" cat "=BEGINS=HERE";
print "$PY_EVAL_SECTION";
print "PY=EVAL=SECTION=ENDS=HERE";

"""

# MAGMA: define the ideal with the non-zero condition. The quotient is the 
# coordinate ring of the Ptolemy variety.
_MAGMA_IDEAL_WITH_NON_ZERO_CONDITION = """
// Setting up the Polynomial ring and ideal

R<$VARIABLES_WITH_NON_ZERO_CONDITION> :=\
 PolynomialRing(RationalField(), $VARIABLE_WITH_NON_ZERO_CONDITION_NUMBER);
MyIdeal := ideal<R |
          $EQUATIONS_WITH_NON_ZERO_CONDITION>;

"""

# MAGMA: call this before computing the decomposition

_MAGMA_PREPARE_DECOMPOSITION = """

// Initialize Q to -1 so that we can check whether an error happend
// by checking that Q is still of type integer.
Q := -1;

// Remember start time to calculate computation time
primTime := Cputime();


"""

# MAGMA: call this after computing and storing the decomposition
# in Q. It will print it.

_MAGMA_PRINT_DECOMPOSITION = """

print "IDEAL=DECOMPOSITION" cat "=TIME: ", Cputime(primTime);

if Type(Q) eq RngIntElt then
    // Some error occured
    print "IDEAL=DECOMPOSITION" cat "=FAILED";
    exit;
else
    // Success
    print "IDEAL=DECOMPOSITION" cat "=BEGINS=HERE";
    Q;
    print "IDEAL=DECOMPOSITION" cat "=ENDS=HERE";


    print "FREE=VARIABLES=IN=COMPONENTS" cat "=BEGINS=HERE";
    N := Names(R);
    isFirstComp := true;
    freeVarStr := "[";
    for Comp in Q do
    
        if isFirstComp then
            isFirstComp := false;
        else
            freeVarStr := freeVarStr cat ",";
        end if;
        freeVarStr := freeVarStr cat "\\n    [ ";
        
        D, Vars := Dimension(Comp);

        isFirstVar := true; 
        for Var in Vars do
            if isFirstVar then
                isFirstVar := false;
            else
                freeVarStr := freeVarStr cat ", ";
            end if;
 
            freeVarStr := freeVarStr cat "\\"" cat N[Var] cat "\\"";
        end for;

        freeVarStr := freeVarStr cat " ]";
    end for;
    freeVarStr := freeVarStr cat "\\n]";
    print freeVarStr;
    print "FREE=VARIABLES=IN=COMPONENTS" cat "=ENDS=HERE";
end if;

print "CPUTIME: ", Cputime(primTime);

"""

########################################
# MAGMA templates

# MAGMA: compute the radical decomposition
# The problem with this magma command is that it sometimes reports the
# dimension just as ">0"

MAGMA_RADICAL_DECOMPOSITION_TEMPLATE = (
    _MAGMA_IDEAL_WITH_NON_ZERO_CONDITION +
    _MAGMA_PRINT_ADDITIONAL_DATA +
    _MAGMA_PREPARE_DECOMPOSITION + """

print "DECOMPOSITION=TYPE: RadicalDecomposition";

Q := RadicalDecomposition(MyIdeal);

""" + _MAGMA_PRINT_DECOMPOSITION)

# MAGMA: compute the primary decomposition of the radical ideal
# This seems to work in the cases where RadicalDecomposition gave
# ">0" as dimension

MAGMA_PRIMARY_DECOMPOSITION_OF_RADICAL_TEMPLATE = (
    _MAGMA_IDEAL_WITH_NON_ZERO_CONDITION +
    _MAGMA_PRINT_ADDITIONAL_DATA +
    _MAGMA_PREPARE_DECOMPOSITION + """

print "DECOMPOSITION=TYPE: Primary Decomposition of Radical";

P, Q := PrimaryDecomposition(Radical(MyIdeal));

""" + _MAGMA_PRINT_DECOMPOSITION)

# MAGMA: give the radicals of the primary ideals in the
# primary decomposition.

# This is different from returning the primary decomposition
# of the radical. Example: I = <x^2, xy>.
# PrimaryDecomposition(I)             = [<x>,<x,y^2>]
# Radicals of PrimaryDecomposition(I) = [<x>,<x,y>]
# Radical(I) = <x>
# PrimaryDecomposition(Radical(I))    = [<x>]
#
# Geometrically, the variety is the y-axis, as a root each point
# on the y-axis has multiplicity 1 except for the origin that has
# multiplicity 2.
#
# For the Ptolemy variety, the implication is this: there is a 0-dimensional
# component in the primary decomposition but it belongs to a
# boundary-unipotent representation that is not boundary-non-degenerate and
# the corresponding cocycle can be deformed.

MAGMA_RADICALS_OF_PRIMARY_DECOMPOSITION_TEMPLATE = (
    _MAGMA_IDEAL_WITH_NON_ZERO_CONDITION +
    _MAGMA_PRINT_ADDITIONAL_DATA +
    _MAGMA_PREPARE_DECOMPOSITION + """

print "DECOMPOSITION=TYPE: Radicals of Primary Decomposition";

P, Q := PrimaryDecomposition(MyIdeal);

""" + _MAGMA_PRINT_DECOMPOSITION)

# MAGMA: the default templated. Our engine to find a solution from the
# Groebner basis works only for prime ideals. See comment of
# MAGMA_RADICALS_OF_PRIMARY_DECOMPOSITION_TEMPLATE why we do not use it.

MAGMA_DEFAULT_TEMPLATE = MAGMA_PRIMARY_DECOMPOSITION_OF_RADICAL_TEMPLATE

########################################
# MAGMA templates for experimentation

# MAGMA: Just compute the Groebner basis.

MAGMA_GROEBNER_BASIS_TEMPLATE = """
R<$VARIABLES> := PolynomialRing(RationalField(), $VARIABLE_NUMBER);
MyIdeal := ideal<R |
          $EQUATIONS>;

""" + _MAGMA_PRINT_ADDITIONAL_DATA + """


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

# MAGMA: use the Variety command

MAGMA_VARIETY_TEMPLATE = """

// Setting up the Polynomial ring and ideal

R<$VARIABLES_WITH_NON_ZERO_CONDITION> :=\
 PolynomialRing(RationalField(), $VARIABLE_WITH_NON_ZERO_CONDITION_NUMBER);
MyIdeal := ideal<R |
          $EQUATIONS_WITH_NON_ZERO_CONDITION>;

""" + _MAGMA_PRINT_ADDITIONAL_DATA + """

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
        print "CPUTIME: ", Cputime(totTime);
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
        print "CPUTIME: ", Cputime(totTime);
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
        print "CPUTIME: ", Cputime(totTime);
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

print "CPUTIME: ", Cputime(totTime);

"""

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
