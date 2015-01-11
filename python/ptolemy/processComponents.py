from .component import NonZeroDimensionalComponent
from . import processFileBase
from . import processRurFile
from . import utilities
from . import coordinates
from .polynomial import Polynomial
from .ptolemyVarietyPrimeIdealGroebnerBasis import PtolemyVarietyPrimeIdealGroebnerBasis

def contains_ideal_components(text):
    return "IDEAL=COMPONENTS=BEGIN" in text

def decomposition_from_components(text):
    
    
    py_eval = processFileBase.get_py_eval(text)
    manifold_thunk = processFileBase.get_manifold_thunk(text)

    variables_section = processFileBase.find_unique_section(
        text, "VARIABLES")

    variables = [ remove_optional_quotes(v.strip())
                  for v in variables_section.split(',') if v.strip() ]

    decomposition = processFileBase.find_unique_section(
        text, "IDEAL=COMPONENTS")

    params, body = processFileBase.extract_parameters_and_body_from_section(
        decomposition)
    
    if not "TYPE" in params.keys():
        raise Exception("No TYPE given for IDEAL=COMPONENTS")

    type = params["TYPE"].strip()

    if not type == "PTOLEMY":
        raise Exception("TYPE '%s' not supported" % type)

    components = processFileBase.find_section(text, "COMPONENT")

    return utilities.MethodMappingList(
        [ process_component(py_eval, manifold_thunk, variables,
                            component)
          for component in components ])

def remove_optional_quotes(s):
    if s[0] in ['"', "'"]:
        assert s[0] == s[-1]
        return s[1:-1]
    return s
    
def process_component(py_eval, manifold_thunk, variables,
                      component):

    params, body = processFileBase.extract_parameters_and_body_from_section(
        component)

    if not "DIMENSION" in params.keys():
        raise Exception("No DIMENSION for COMPONENT of IDEAL=COMPONENTS")

    dimension = int(params["DIMENSION"].strip())

    free_variables = None

    if dimension > 0 and "FREE=VARIABLES" in params.keys():
        free_vars_str = processFileBase.remove_optional_outer_square_brackets(
            params["FREE=VARIABLES"])

        free_variables = [ remove_optional_quotes(v.strip())
                           for v in free_vars_str.split() ]

    if dimension == 0:
        return process_solutions_provider(py_eval, manifold_thunk, body, 0,
                                          variables)
    else:
        witnesses = []
        witnesses_sections = processFileBase.find_section(body, "WITNESSES")
        if witnesses_sections:
            assert len(witnesses_sections) == 1
            witnesses = process_witnesses(
                py_eval, manifold_thunk, witnesses_sections[0], dimension,
                variables)
        
        # process Groebner Basis ???, add to NonZeroDimensionalComponent ???

        return NonZeroDimensionalComponent(
            witnesses = witnesses,
            dimension = dimension,
            free_variables = free_variables)

def process_witnesses(py_eval, manifold_thunk, witnesses_section, for_dimension,
                      variables):
    
    return [ process_solutions_provider(
            py_eval, manifold_thunk, witness_section, for_dimension, variables)
             for witness_section
             in processFileBase.find_section(witnesses_section, "WITNESS") ]

def process_solutions_provider(py_eval, manifold_thunk, text, for_dimension,
                               variables):

    rur_section = processFileBase.find_section(text, "MAPLE=LIKE=RUR")
    if rur_section:
        assert len(rur_section) == 1

        return SolutionContainer(
            coordinates.PtolemyCoordinates(
                processRurFile.parse_maple_like_rur(rur_section[0]),
                is_numerical = False,
                py_eval_section = py_eval,
                manifold_thunk = manifold_thunk))

    gb_section = processFileBase.find_section(text, "GROEBNER=BASIS")

    if gb_section:
        assert len(gb_section) == 1

        params, body = (
            processFileBase.extract_parameters_and_body_from_section(
                gb_section[0]))

        body = processFileBase.remove_optional_outer_square_brackets(
            utilities.join_long_lines_deleting_whitespace(body))
        
        polys = [ Polynomial.parse_string(p)
                  for p in body.replace('\n', ' ').split(',') ]

        if not 'TERM=ORDER' in params.keys():
            raise Exception("No term order given for Groebner basis")

        term_order = params["TERM=ORDER"].strip().lower()
        
        return PtolemyVarietyPrimeIdealGroebnerBasis(
            polys = polys,
            term_order = term_order,
            size = None,
            dimension = 0,
            is_prime = True,
            free_variables = None,
            py_eval = py_eval,
            manifold_thunk = manifold_thunk)

    rs_rur = processFileBase.find_section(text, "RS=RUR")

    if rs_rur:
        assert len(rs_rur) == 1

        return SolutionContainer(
            coordinates.PtolemyCoordinates(
                processRurFile.parse_rs_rur(rs_rur[0], variables),
                is_numerical = False,
                py_eval_section = py_eval,
                manifold_thunk = manifold_thunk))
    
    raise Exception("No parsable solution type given: %s..." % text[:100])

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
