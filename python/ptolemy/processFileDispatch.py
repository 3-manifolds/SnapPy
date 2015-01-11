"""
This module takes a textual representation of the solutions to a
Ptolemy variety and parses it by detecting the type of the textural
representation and dispatching it to the corresponding module.
"""

from . import processMagmaFile
from . import processRurFile
from . import processComponents

def parse_decomposition(text):

    if processMagmaFile.contains_magma_output(text):
        return processMagmaFile.decomposition_from_magma(text)

    # Deprecated, will be removed soon!
    if processRurFile.contains_rur(text):
        return processRurFile.decomposition_from_rur(text)

    if processComponents.contains_ideal_components(text):
        return processComponents.decomposition_from_components(text)

    raise Exception("Solution file format not recognized")
    
def parse_decomposition_from_file(filename):
    
    return parse_decomposition(open(filename).read())

def parse_solutions(text, numerical = False):

    """
    Reads the text containing the solutions from a magma computation
    or a rur computation and returns a list of solutions.
    A non-zero dimensional component of the variety is reported as
    NonZeroDimensionalComponent.
    """

    return parse_decomposition(text).solutions(numerical)

    
def parse_solutions_from_file(filename, numerical = False):
    
    """
    As parse_solutions, but takes a filename instead.
    """

    return parse_solutions(open(filename).read(), numerical)
