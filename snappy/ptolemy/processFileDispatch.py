"""
This module takes a textual representation of the solutions to a
Ptolemy variety and parses it by detecting the type of the textural
representation and dispatching it to the corresponding module.
"""

import processMagmaFile

def get_solutions(text, numerical = False):

    """
    Reads the text containing the solutions from a magma computation
    and returns a list of solutions.
    A non-zero dimensional component of the variety is reported as
    NonZeroDimensionalComponent.
    """

    if processMagmaFile.contains_magma_output(text):
        return processMagmaFile.solutions_from_magma(text, numerical)

    raise Exception("Solutions format not recognized")

    
def get_solutions_from_file(filename, numerical = False):
    
    """
    As solutions_from_text, but takes a filename instead.
    """

    return get_solutions(open(filename).read(), numerical)
