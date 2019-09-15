from snappy.verify.complexVolume.extended_bloch import (
    compute_complex_volume_from_lifted_ptolemys)

from snappy.dev import extended
from snappy.dev import giac_rur

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
    a list of list of dictionaries assinging complex intervals to
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

def adjust_meridian_and_longitude(M, m, l, precision = 53):
    """
    Given a closed manifold and any log of the holonomy of the meridian and
    longitude, adjust logs by multiplies of 2 pi i such that the peripheral
    curves goes to 0.
    """

    RIF = RealIntervalField(precision)
    two_pi = RIF(2*pi)

    # (m_fill, l_fill) Dehn-filling
    m_fill, l_fill = [int(x) for x in M.cusp_info()[0]['filling']]

    # Compute what the peripheral curves goes to right now
    p_interval = (m_fill * m + l_fill * l).imag() / two_pi
    is_int, p = p_interval.is_int()

    if not is_int:
        raise Exception("Expected multiple of 2 * pi * i (increase precision?)")

    if p == 0:
        # Nothing to do
        return m, l

    # Compute by what multiple of 2 pi i to adjust
    g, a, b = xgcd(m_fill, l_fill)
    m -= p * a * two_pi * sage.all.I
    l -= p * b * two_pi * sage.all.I

    # For sanity, double check that we compute it right.
    p_interval = (m_fill * m + l_fill * l).imag() / two_pi
    is_int, p = p_interval.is_int()

    if not is_int:
        raise Exception("Expected multiple of 2 * pi * i (increase precision?)")

    if p != 0:
        # Nothing to do
        raise Exception("Expected 0")

    return m, l

def lift_ptolemy_coordinates(M, solution, full_var_dict, precision = 53):
    """
    Given a closed manifold (as Dehn-filling on 1-cusped manifold) and an
    assignment of subset of ptolemy variables and the full var dict, compute
    logs for all Ptolemy's.
    """

    lifted = { str(k) : log(v)
               for k, v in solution.items()
               if str(k)[0].islower() }

    m, l = adjust_meridian_and_longitude(
        M, lifted['m'], lifted['l'], precision)

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
        [ compute_complex_volume_from_lifted_ptolemys(
                M.num_tetrahedra(),
                lift_ptolemy_coordinates(M, sol, full_var_dict, precision))
          for sol in galois_conjugates ]
        for galois_conjugates in representative_ptolemys ]

def has_value(v, values):

    RIF = RealIntervalField(212)

    for value in values:
        if abs(RIF(v.imag()) - RIF(value.imag())) < 1e-20:
            r = (RIF(v.real()) - RIF(value.real())) / RIF(pi**2/6)

            is_int, k = r.is_int()
            if is_int:
                if abs(r - k) < 1e-20:
                    return True
    return False

if __name__ == '__main__':

    from snappy import ManifoldHP

    M = ManifoldHP("5_2")
    M.chern_simons()
    M.dehn_fill((1,2))

    cvol = M.complex_volume() * sage.all.I

    # Because of ordering issues, only correct up to pi^2/6
    cvols = sum(complex_volumes(M, precision = 300), [])

    print(has_value(cvol, cvols))
