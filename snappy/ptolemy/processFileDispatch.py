"""
This module takes a textual representation of the solutions to a
Ptolemy variety and parses it by detecting the type of the textural
representation and dispatching it to the corresponding module.
"""

from . import processMagmaFile
from . import processRurFile

def parse_solutions(text, numerical = False):

    """
    Reads the text containing the solutions from a magma computation
    or a rur computation and returns a list of solutions.
    A non-zero dimensional component of the variety is reported as
    NonZeroDimensionalComponent.
    """

    if processMagmaFile.contains_magma_output(text):
        return processMagmaFile.solutions_from_magma(text, numerical)

    if processRurFile.contains_rur(text):
        return processRurFile.solutions_from_rur(text, numerical)

    raise Exception("Solutions format not recognized")

    
def parse_solutions_from_file(filename, numerical = False):
    
    """
    As parse_solutions, but takes a filename instead.
    """

    return parse_solutions(open(filename).read(), numerical)
