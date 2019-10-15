from ...sage_helper import _within_sage, sage_method

if _within_sage:
    from sage.all import (ComplexBallField,
                          RealField,
                          Integer, exp, pi)

    import sage.all

@sage_method
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

@sage_method
def compute_p_from_w_and_parity(w, parity):
    """
    Compute p such that w - p * pi * i should have imaginary part between
    -pi and pi and p has the same parity as the given value for parity
    (the given value is supposed to be 0 or 1).

    Note that this computation is not verified.
    """

    RF = RealField(w.parent().precision())
    real_part = (w.imag().center() / RF(pi) - parity) / 2
    return 2 * Integer(real_part.round()) + parity

@sage_method
def compute_z_p_q_from_flattening_w0_w1(w0, w1):
    """
    Given w0 and w1 such that +- exp(w0) +- exp(-w1) = 1, compute
    a triple [z; p, q] such that
    w0 = log(z) + p * pi * i and w1 = -log(1-z) + q * pi * i.
    
    While z is and the parities of p and q are verified, p and q are
    not verified in the following sense:
    w0 - p * pi * i and w1 + q * pi * i are likely to have imaginary
    part between -pi and pi, but this is not verified.
    """

    z, p_parity, q_parity = compute_z_and_parities_from_flattening_w0_w1(w0, w1)

    return (z,
            compute_p_from_w_and_parity(w0, p_parity),
            compute_p_from_w_and_parity(w1, q_parity))

@sage_method
def my_dilog(z):
    """
    Compute dilogarithm using complex ball field.
    The dilogarithm isn't implemented for ComplexIntervalField itself, so
    we use ComplexBallField. Note that ComplexBallField is conservative
    about branch cuts. For Li_2(2+-i * epsilon), it returns the interval
    containing both Li_2(2+i * epsilon) and Li_2(2-i * epsilon).

    Thus, we need to avoid calling this function with a value near real numbers
    greater 1.
    """

    CIF = z.parent()
    CBF = ComplexBallField(CIF.precision())

    return CIF(CBF(z).polylog(2))

@sage_method
def is_imaginary_part_bounded(z, v):
    """
    Check that the imaginary part of z is in (-v, v).
    """

    imag = z.imag()
    return -v < imag and imag < v

@sage_method
def compute_Neumanns_Rogers_dilog_from_flattening_w0_w1(w0, w1):
    """
    Given a flattening w0, w1 such that +- exp(w0) +- exp(-w1) = 1, compute
    the complex volume given by R(z;p,q) (equation before Proposition 2.5 in
    Neumann's Extended Bloch group and the Cheeger-Chern-Simons class).
    """

    RIF = w0.parent().real_field()
    my_pi = RIF(pi)

    # Compute [z; p, q]
    z, p, q = compute_z_p_q_from_flattening_w0_w1(w0, w1)

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
    # Note that we should always be able to avoid this by increasing bits_prec
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

        return ( term1 + term2) / 2 + my_dilog(z) - my_pi * my_pi / 6
    else:
        # Check that we can apply equation (2)
        if not is_imaginary_part_bounded(logZ, 2):
            raise Exception("Problem with computig Neumanns dilog using (2)",
                            z, logZ)

        return (-term1 + term2) / 2 - my_dilog(1 - z)

@sage_method
def compute_complex_volume_of_simplex_from_lifted_ptolemys(index, ptolemys):
    """
    Given lifted Ptolemy coordinates for a triangulation (as dictionary),
    compute the complex volume contribution by the simplex with given index.
    """

    # The six Ptolemy coordinates for the given simplex
    c_1100 = ptolemys['c_1100_%d' % index]
    c_1010 = ptolemys['c_1010_%d' % index]
    c_1001 = ptolemys['c_1001_%d' % index]
    c_0110 = ptolemys['c_0110_%d' % index]
    c_0101 = ptolemys['c_0101_%d' % index]
    c_0011 = ptolemys['c_0011_%d' % index]

    # Compute Neumann's flattening (w0, w1) from Ptolemy coordinates
    w0 = c_1010 + c_0101 - c_1001 - c_0110
    w1 = c_1001 + c_0110 - c_1100 - c_0011

    # Compute Neumann's version of Roger's dilogarithm from flattening.
    return compute_Neumanns_Rogers_dilog_from_flattening_w0_w1(w0, w1)

@sage_method
def compute_complex_volume_from_lifted_ptolemys_no_torsion_adjustment(
        num_tetrahedra, ptolemys):
    """
    Given lifted Ptolemy coordinates for a triangulation (as dictionary)
    and the number of tetrahedra, compute the complex volume (where
    the real part is the Chern-Simons and the imaginary part is the
    volume).

    This sums of the dilogs across tetrahedra without adjusting for the
    fact that the triangulation might not be ordered.
    Thus, the Chern-Simons is correct only up to multiples of pi^2/6.
    """

    return sum(
        [ compute_complex_volume_of_simplex_from_lifted_ptolemys(
                index, ptolemys)
          for index in range(num_tetrahedra) ])
