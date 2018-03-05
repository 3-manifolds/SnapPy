import re

from .utilities import MethodMappingList
from .component import NonZeroDimensionalComponent
from .coordinates import PtolemyCoordinates
from .rur import RUR
from . import processFileBase
from ..pari import pari

def contains_rur(text):
    return 'RUR=DECOMPOSITION=BEGINS' in text

def decomposition_from_rur(text):

    py_eval = processFileBase.get_py_eval(text)
    manifold_thunk = processFileBase.get_manifold_thunk(text)

    rursection = processFileBase.find_unique_section(text, "RUR=DECOMPOSITION")

    rurs = processFileBase.find_section(rursection, "COMPONENT")

    result = MethodMappingList(
        [ SolutionContainer(
                _process_rur_component(rur, py_eval, manifold_thunk))
          for rur in rurs ])

    return result

def _process_rur_component(text, py_eval, manifold_thunk):
    lines = text.split('\n')

    dimension = None
    is_prime = False
    format = None

    body = ""
    for line in lines:
        m = re.match(r'==(.*):(.*)', line)
        if m:
            key, val = m.groups()
            val = val.strip()
            if key == "DIMENSION":
                dimension = int(val)
            elif key == "IS_PRIME":
                if val == "TRUE":
                    is_prime = True
                elif val != "FALSE":
                    raise Exception("IS_PRIME needs to be TRUE or FALSE")
            elif key == "FORMAT":
                format = val
            else:
                raise Exception("Unrecognized key %s" % key)
        else:
            body += line + "\n"
   
    if dimension is None:
        return NonZeroDimensionalComponent()
    if dimension > 0:
        return NonZeroDimensionalComponent(dimension = dimension)

    if format is None:
        raise Exception("No format specified")

    if format == "MAPLE-LIKE":
        d = parse_maple_like_rur(body.strip())
        return PtolemyCoordinates(d, is_numerical = False,
                                  py_eval_section = py_eval,
                                  manifold_thunk = manifold_thunk)
    else:
        raise Exception("Unknown format %s" % format)


def parse_maple_like_rur(text):

    m = re.match("(.*?)\s*=\s*0\s*,\s*\{(.*?)\}", text, re.DOTALL)

    if not m:
        raise Exception("Format not detected")

    poly_text, assignments_text = m.groups()
    var = _find_var_of_poly(poly_text)

    poly = pari(poly_text.replace(var, 'x'))
    assignments_text = assignments_text.replace(var, 'x')

    return dict([_parse_assignment(assignment, poly)
                 for assignment in assignments_text.split(',')])

def parse_rs_rur(text, variables):
    
    m = re.match(r'\[([^,\]]+),\s*([^,\]]+),\s*\[([^\]]+)\]\s*,\s*\[\s*\]\s*\]',
                 text, re.DOTALL)

    if not m:
        raise Exception("Format not detected")

    extension_str, denominator_str, numerators_str = m.groups()

    var = _find_var_of_poly(extension_str)

    extension = pari(extension_str.replace(var, 'x'))
    denominator = pari(denominator_str.replace(var, 'x'))
    numerators = [ pari(numerator_str.replace(var, 'x'))
                   for numerator_str in numerators_str.split(',') ]
    
    fracs = [
        RUR( [ (  numerator.Mod(extension),  1),
               (denominator.Mod(extension), -1) ] )
        for numerator in numerators ]

    return dict(zip(variables, fracs) + [('1', RUR.from_int(1))])
    
def _find_var_of_poly(text):
    return re.search(r'[_A-Za-z][_A-Za-z0-9]*',text).group(0)

def _parse_assignment(text, poly):
    var, expression = re.match(r'\s*([_A-Za-z][_A-Za-z0-9]*)\s*=\s*(.*)$',
                               text).groups()
    
    return (
        var,
        RUR.from_pari_fraction_and_number_field(pari(expression), poly))

_test_case = """
_Z^6-3*_Z^5+3*_Z^4-2*_Z^3+_Z-1 = 0,
 {c_0012_0 = 1, c_0012_2 = _Z^5-3*_Z^4+3*_Z^3-_Z^2-2*_Z+1, c_0111_2 = -_Z*(_Z^4-2*_Z^3+_Z-1), c_1011_0 = -_Z^2*(-1+_Z), c_1101_0 = -_Z*(-1+_Z), c_1101_1 = _Z*(-1+_Z), c_1110_0 = _Z, a = (_Z^3+_Z^2)/(_Z^3+34*_Z^2+3)}
"""


class SolutionContainer():
    def __init__(self, solutions):
        self._solutions = solutions

    def solutions(self, numerical = False):
        if numerical:
            return self._solutions.numerical()
        else:
            return self._solutions

    def number_field(self):
        return self._solutions.number_field()

# print parse_maple_like_rur(_test_case)
