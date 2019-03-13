from snappy.dev import extended
from snappy.dev import giac_rur

from sage.all import (RealIntervalField, ComplexIntervalField,
                      RealBallField, ComplexBallField,
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
    
    # Should be called univariate
    rur = giac_rur.rational_unimodular_representation(I)

    return [
        evaluate_at_roots(numberField, exact_values, precision)
        for numberField, exact_values in rur ], full_var_dict

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

def compute_z_from_flattening_w0_w1(w0, w1):
    """
    Given a pair (w0, w1) with +- exp(w0) +- exp(-w1) = 1, compute z such that
    z = +-exp(w0) and 1/(1-z) = +- exp(w1).
    """

    e0 = exp( w0)
    e1 = exp(-w1)
    
    l = [ sign0 * e0
          for sign0 in [+1, -1]
          for sign1 in [+1, -1]
          if Integer(1) in sign0 * e0 + sign1 * e1 ]
    if not len(l) == 1:
        raise Exception("Bad flattening %s %s %s" % (w0, w1, len(l)))

    return l[0]

def my_dilog(z, precision = 53):
    """
    Compute dilogarithm using complex ball field.
    It is unfortunate that dilogarithm isn't implemented for
    ComplexIntervalField itself. In particular, since ComplexBallField returns
    a large interval when near a branch cut instead of lifting.
    """

    CBF = ComplexBallField(precision)
    CIF = ComplexIntervalField(precision)
    
    return CIF(CBF(z).polylog(2))

def compute_Neumanns_Rogers_dilog_from_flattening_w0_w1(w0, w1, precision = 53):
    """
    Given a flattening w0, w1 such that +- exp(w0) +- exp(-w1) = 1, compute
    the complex volume computed by L(z;p,q).
    """

    RIF = RealIntervalField(precision)
    my_pi = RIF(pi)

    z = compute_z_from_flattening_w0_w1(w0, w1)
    
    logZ = log(z)
    logOneMinusZ = log(1 - z)

    p_interval = (w0 - logZ).imag() / my_pi
    is_int, p = p_interval.is_int()
    if not is_int:
        raise Exception("Expected integer for p")

    q_interval = (w1 + logOneMinusZ).imag() / my_pi
    is_int, q = q_interval.is_int()
    if not is_int:
        raise Exception("Expected integer for q")

    t1 = logZ + p * my_pi * sage.all.I
    t2 = logOneMinusZ + q * my_pi * sage.all.I

    return my_dilog(z, precision) + t1 * t2 / 2 - my_pi * my_pi / 6

def compute_complex_volume_of_simplex_from_lifted_ptolemys(index, ptolemys,
                                                           precision):
    c_1100 = ptolemys['c_1100_%d' % index]
    c_1010 = ptolemys['c_1010_%d' % index]
    c_1001 = ptolemys['c_1001_%d' % index]
    c_0110 = ptolemys['c_0110_%d' % index]
    c_0101 = ptolemys['c_0101_%d' % index]
    c_0011 = ptolemys['c_0011_%d' % index]
    
    w0 = c_1010 + c_0101 - c_1001 - c_0110
    w1 = c_1001 + c_0110 - c_1100 - c_0011
    
    return compute_Neumanns_Rogers_dilog_from_flattening_w0_w1(
        w0, w1, precision)

def compute_complex_volume_from_lifted_ptolemys(M, ptolemys, precision):
    return sum(
        [ compute_complex_volume_of_simplex_from_lifted_ptolemys(
                index, ptolemys, precision)
          for index in range(M.num_tetrahedra()) ])
    

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
                M,
                lift_ptolemy_coordinates(M, sol, full_var_dict, precision),
                precision)
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

    print has_value(cvol, cvols)
    
