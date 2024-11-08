from snappy.verify.complex_volume.adjust_torsion import (
    verified_complex_volume_from_lifted_ptolemys)
from snappy.verify.complex_volume.closed import zero_lifted_holonomy

from snappy.dev.extended_ptolemy import extended
from snappy.dev.extended_ptolemy import giac_rur

import snappy.snap.t3mlite as t3m

from sage.all import (RealIntervalField, ComplexIntervalField,
                      RealBallField, ComplexBallField,
                      RealField,
                      Integer,
                      prod, log, exp, pi, xgcd)
import sage.all

import re

def evaluate_at_roots(numberField, exact_values, precision = 53):
    """
    numberField is a sage number field.
    exact_values a dictionary where values are elements in that number field.
    precision is desired precision in bits.

    For each embedding of the number field, evaluates the values in the
    dictionary and produces a new dictionary with the same keys.
    The new dictionaries are returned in a list.
    """

    CIF = ComplexIntervalField(precision)
    return [{ k : v.lift().substitute(root) for k, v in exact_values.items()}
            for root, multiplicity in numberField.polynomial().roots(CIF)]

def compute_representative_ptolemys_and_full_var_dict(M, precision = 53):
    """
    Given a closed manifold (as Dehn-filling on 1-cusped manifold), compute
    a list of list of dictionaries assigning complex intervals to
    a subset of ptolemy variables and the full var dictionary to expand these
    to all variables.
    (Outer list for components, inner list for Galois conjugates).
    """

    I, full_var_dict = extended.ptolemy_ideal_for_filled(
        M, return_full_var_dict = 'data', notation = 'full')

    rur = giac_rur.rational_univariate_representation(I)

    return [
        evaluate_at_roots(numberField, exact_values, precision)
        for numberField, exact_values, mult in rur ], full_var_dict

def lift_ptolemy_coordinates(M, solution, full_var_dict):
    """
    Given a closed manifold (as Dehn-filling on 1-cusped manifold) and an
    assignment of subset of ptolemy variables and the full var dict, compute
    logs for all Ptolemy's.
    """

    lifted = { str(k) : log(v)
               for k, v in solution.items()
               if str(k)[0].islower() }

    m, l = zero_lifted_holonomy(
        M, lifted['m'], lifted['l'], 2)

    return { k : lifted[name] - (m_count * m + l_count * l)
             for k, (sign, m_count, l_count, name) in full_var_dict.items()
             if k[0] == 'c'}

def complex_volumes(M, precision = 53):
    """
    Compute all complex volumes from the extended Ptolemy variety for the
    closed manifold M (given as Dehn-filling on 1-cusped manifold).
    Note: not every volume might correspond to a representation factoring
    through the closed manifold. In particular, we get the complex volume
    of the geometric representation of the cusped manifold.
    """

    representative_ptolemys, full_var_dict = (
        compute_representative_ptolemys_and_full_var_dict(M, precision))

    return [
        [ verified_complex_volume_from_lifted_ptolemys(
                t3m.Mcomplex(M),
                lift_ptolemy_coordinates(M, sol, full_var_dict))
          for sol in galois_conjugates ]
        for galois_conjugates in representative_ptolemys ]

def has_value(v, values):

    RIF = RealIntervalField(212)

    for value in values:
        if abs(RIF(v.imag()) - RIF(value.imag())) < RIF(1e-20):
            r = (RIF(v.real()) - RIF(value.real())) / RIF(pi**2/2)

            is_int, k = r.is_int()
            if is_int:
                if abs(r - k) < RIF(1e-20):
                    return True
    return False

if __name__ == '__main__':
    # Use 5_2(1,2) as an example.
    from snappy import ManifoldHP

    M = ManifoldHP("5_2")
    M.chern_simons()
    M.dehn_fill((1,2))

    # Complex volume of geometric representation of 5_2(1,2).
    # Unverified.
    geom_cvol = M.complex_volume()

    # This file uses a different convention (simply adding Neumann's
    # extensions of Roger's dilogarithm to the Abelian cover of C \ {0,1})
    # from the rest of SnapPy, so convert.
    geom_cvol = sage.all.I * geom_cvol

    # The complex volumes as verified intervals of the SL(2,C)-representations
    # detected by the (spun-)triangulation.
    # Grouped by algebraic component of the extended Ptolemy variety.
    #
    # We can only compute it in this file up to multiples of pi^2/2.
    #
    # (Note that Neumann's method can compute it up to multiple of pi^2 but
    # assumes an ordered triangulation. Most census triangulations are not
    # orderable and would require subdivision).
    #
    # The method uses SageMath's giac to compute the rational univariate
    # representation and evaluates it at the roots of the defining polynomial
    # using complex intervals.
    #
    # It can be accessed from sage (with SnapPy installed) by:
    # from snappy.dev.extended_ptolemy.complexVolumesClosed import complex_volumes
    #
    # Tested with SageMath 9.7 and 10.0.
    #
    cvols_by_component = complex_volumes(M, precision = 300)

    print(cvols_by_component)

    # Concatenate the grouped complex volumes to just have a single list.
    cvols = sum(cvols_by_component, [])

    if not has_value(geom_cvol, cvols):
        raise RuntimeError(
            "There is a problem. The complex volume of the geometric "
            "representation is not included.")
