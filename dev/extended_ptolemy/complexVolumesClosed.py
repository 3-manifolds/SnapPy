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

def compute_z_and_parities_from_flattening_w0_w1(w0, w1):
    """
    Given a pair (w0, w1) with +- exp(w0) +- exp(-w1) = 1, compute (z, p, q)
    such that z = (-1)^p * exp(w0) and 1/(1-z) = (-1)^q exp(w1)
    where p, q in {0,1}.
    """

    e0 = exp( w0)
    e1 = exp(-w1)

    l = [ (((-1) ** p) * e0, p, q)
          for p in [ 0, 1]
          for q in [ 0, 1]
          if Integer(1) in ((-1) ** p) * e0 + ((-1) ** q) * e1 ]
    if not len(l) == 1:
        raise Exception("Bad flattening %s %s %s" % (w0, w1, len(l)))

    return l[0]

def compute_p_from_w_and_parity(w, parity, precision):
    """
    Compute p such that w - p * pi * i should have imaginary part between
    -pi and pi and p has the same parity as the given value for parity
    (the given value is supposed to be 0 or 1).

    Note that this computation is not verified.
    """

    RF = RealField(precision)
    real_part = (w.imag().center() / RF(pi) - parity) / 2
    return 2 * Integer(real_part.round()) + parity

def compute_z_p_q_from_flattening_w0_w1(w0, w1, precision):
    """
    Given w0 and w1 such that +- exp(w0) +- exp(-w1) = 1, compute
    a triple [z; p, q] such that
    w0 = log(z) + p * pi * i and w1 = -log(1-z) + q * pi * i.
    
    While z is and the parities of p and q are verified, p and q are
    not verified in the following sense:
    w0 - p * pi * i and w1 + q * pi * i are likely to have imaginary
    part between -pi and pi, but this is not verified.
    """

    RF = RealField(precision)
    my_pi = RF(pi)

    z, p_parity, q_parity = compute_z_and_parities_from_flattening_w0_w1(w0, w1)

    return (z,
            compute_p_from_w_and_parity(w0, p_parity, precision),
            compute_p_from_w_and_parity(w1, q_parity, precision))

def my_dilog(z, precision = 53):
    """
    Compute dilogarithm using complex ball field.
    The dilogarithm isn't implemented for ComplexIntervalField itself, so
    we use ComplexBallField. Note that ComplexBallField is conservative
    about branch cuts. For Li_2(2+-i * epsilon), it returns the interval
    containing both Li_2(2+i * epsilon) and Li_2(2-i * epsilon).

    Thus, we need to avoid calling this function with a value near real numbers
    greater 1.
    """

    CBF = ComplexBallField(precision)
    CIF = ComplexIntervalField(precision)

    return CIF(CBF(z).polylog(2))

def is_imaginary_part_bounded(z, v):
    """
    Check that the imaginary part of z is in (-v, v).
    """

    imag = z.imag()
    return -v < imag and imag < v

def compute_Neumanns_Rogers_dilog_from_flattening_w0_w1(w0, w1, precision = 53):
    """
    Given a flattening w0, w1 such that +- exp(w0) +- exp(-w1) = 1, compute
    the complex volume given by R(z;p,q) (equation before Proposition 2.5 in
    Neumann's Extended Bloch group and the Cheeger–Chern–Simons class).
    """

    RIF = RealIntervalField(precision)
    my_pi = RIF(pi)

    # Compute [z; p, q]
    z, p, q = compute_z_p_q_from_flattening_w0_w1(w0, w1, precision)

    # Note that the values computed for log(z) and log(1-z)
    # are not verified to have the imaginary part between -pi and pi.
    logZ         =    w0 - my_pi * p * sage.all.I
    logOneMinusZ = - (w1 - my_pi * q * sage.all.I)
    
    # Neumann's formula for the complex volume is
    #
    # (1) R(z; p, q) =   Li_2(  z) + ( term1 + term2) / 2 - pi^2/6
    #
    # where
    #     term1 = log(z) * log(1-z)
    #     term2 = pi * i * (p * log(1-z) + q * log(z))
    #
    # Using Li_2(z) + Li_1(1-z) = pi^2/6 - log(z) * log(1-z), we also get
    #
    # (2) R(z; p, q) = - Li_2(1-z) + (-term1 + term2) / 2
    #
    # We use (1) when Re(z) < 1/2 and (2) otherwise.
    #
    # Note that if we use (1), we do not rely on the value computed for log(z)
    # to have imaginary part between -pi and pi (because p was not computed
    # such that we have this property). More precisely, if we add 2 to p, the
    # value computed for log(z) changes by -2 * pi * i, so term1 changes by
    # -2 * pi * i * log(1-z) but this is compensated by the change in term2.
    # We do, however, need to check that the value of log(1-z) has imaginary
    # part between -pi and pi. Since We have Re(z) < 1/2, we indeed expect that
    # the imaginary part is between -pi/2 and pi/2 and can check the stronger
    # condition that the imaginary part is between -2 and 2.
    # We need to make sure that Li_2(z) is evaluated correctly. We always
    # want to take the main branch. If z is close to the branch cut ([1,inf)),
    # the choice is ambiguous but we are safe since my_dilog would
    # conservatively return the large interval containing both branch choices.
    # Note that we should always be able to avoid this by increasing precision
    # since we only use (1) if the interval for z is centered to the left of
    # the line with real part 1/2.
    #
    # Similar considerations apply to (2) used when Re(z) > 1/2.

    term1 = logZ * logOneMinusZ
    term2 = my_pi * sage.all.I * (p * logOneMinusZ + q * logZ)

    if z.real().center() < 0.5:
        # Check that we can apply equation (1)
        if not is_imaginary_part_bounded(logOneMinusZ, 2):
            raise Exception("Problem with computig Neumanns dilog using (1)",
                            z, logOneMinusZ)

        return ( term1 + term2) / 2 + my_dilog(z, precision) - my_pi * my_pi / 6
    else:
        # Check that we can apply equation (2)
        if not is_imaginary_part_bounded(logZ, 2):
            raise Exception("Problem with computig Neumanns dilog using (2)",
                            z, logZ)

        return (-term1 + term2) / 2 - my_dilog(1 - z, precision)

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

    print(has_value(cvol, cvols))
