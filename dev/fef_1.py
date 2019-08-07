from snappy import *

def _coprime(a, b):
    if b == 0:
        return a == 1 or a == -1
    return _coprime(b, a % b)

def slopes_of_length_leq_6_for_cusp(t):
    m_trans, l_trans = t

    RF = l_trans.parent()

    return _interval_slopes_of_length_leq_6_from_translations(t)

def _max_int_in_interval(i):
    return i.upper().floor()

def _min_int_in_interval(i):
    return i.lower().ceil()

def _interval_slopes_of_length_leq_6_from_translations(t):
    m_trans, l_trans = t

    RIF = l_trans.parent()

    six = RIF(6)
    max_abs_m = _max_int_in_interval(six / abs(m_trans.imag()))

    for m in range(0, max_abs_m + 1):
        total_m_trans = m * m_trans
 
        max_real_range_sqr = (six ** 2 - total_m_trans.imag() ** 2).upper()

        if max_real_range_sqr > 0:
            max_real_range = RIF(max_real_range_sqr).sqrt()

            min_l = _min_int_in_interval(
                (- total_m_trans.real() - max_real_range) / l_trans.real())
            max_l = _max_int_in_interval(
                (- total_m_trans.real() + max_real_range) / l_trans.real())

            if m == 0:
                min_l = 1

            for l in range(min_l, max_l + 1):
                if _coprime(m, l):
                    yield (m, l)

M = Manifold("m004")

a = M.cusp_translations(verified=True)[0]

for m, l in slopes_of_length_leq_6_from_translations(a):
    print m, l
