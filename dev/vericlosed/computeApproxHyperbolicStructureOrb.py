from .parseVertexGramMatrixFile import *
from .verificationError import *

from .orb import __path__ as orb_path

from snappy.snap.t3mlite import Mcomplex

import subprocess
import tempfile
import shutil
import os
    
__all__ = ['compute_approx_hyperbolic_structure_orb']

class TmpDir(object):
    """
    To be used in a with statement, creating a temporary
    directory deleted at the end of the with statement.
    """

    def __init__(self, delete = True):
        self.delete = delete

    def __enter__(self):
        self.path = tempfile.mkdtemp()
        return self

    def __exit__(self, type, value, traceback):
        if self.delete:
            shutil.rmtree(self.path)

def _compute_vertex_gram_matrix_file_orb(triangulation, tmpDir, verbose = False):
    """
    Calls Orb to compute an approximate solution to the edge equation
    for the given snappy.Triangulation (using directory at path tmpDir).
    """

    # Path for triangulation
    basename = os.path.join(tmpDir, 'manifold.tri')

    # Save triangulation to disk
    triangulation.save(basename)

    if verbose:
        print("Calling orb to find approximate solution on", basename)

    # Find Orb wrapper
    path = os.path.join(orb_path[0], 'orb_solution_for_snappea_finite_triangulation_mac')

    try:
        if verbose:
            stdout = None
        else:
            # Redirect output to nirvana
            stdout = open(os.devnull, 'w')

        # Call Orb
        ret = subprocess.call([path, basename], stdout = stdout)
    except OSError as e:
        raise OrbMissingError(
            "Seems orb binary %s is missing. Original exception: %s" % (
                path, e))

    if ret != 0:
        if ret <= 7:
            raise OrbSolutionTypeError("Solution type orb found was", ret)
        raise UnknownOrbFailureError("Exit code", ret)

    # Path of output that Orb produced
    return basename + '.vgm'

def compute_approx_hyperbolic_structure_orb(triangulation, verbose = False):
    """
    Calls Orb to compute an approximate solution to the edge equation
    for the given snappy.Triangulation.

    The result is a veriClosed.HyperbolicStructure where the edge lengths
    are in SageMath's RealDoubleField.
    """

    # Input and output to Orb are in a temporary directory deleted
    # at the end
    with TmpDir() as tmp_dir:
        # Invoke Orb to compute file containing vertex gram matrices
        vgm_file_path = _compute_vertex_gram_matrix_file_orb(
            triangulation, tmp_dir.path, verbose = verbose)
        # Parse the file to obtain hyperbolic structure
        return compute_approx_hyperbolic_structure_from_vertex_gram_matrix_file(
            Mcomplex(triangulation), vgm_file_path)
