from sage.all import sqrt, gcd

def short_slopes_from_cusp_shape_and_area(
            cusp_shape, cusp_area, length = 6):
    """
    cusp_shape is longitude over meridian (conjugate).
    l/m
    """

    return short_slopes_from_translations(
        translations_from_cusp_shape_and_area(cusp_shape, cusp_area),
        length)

def translations_from_cusp_shape_and_area(
            cusp_shape, cusp_area):
    
    scale = sqrt(cusp_area / cusp_shape.imag())
    return (scale, cusp_shape * scale)

def short_slopes_from_translations(translations, length = 6):
    m_tran, l_tran = translations

    RIF = m_tran.parent()

    length = RIF(length)

    result = []
    
    max_abs_l = _max_int_in_interval(length / abs(l_tran.imag()))

    for l in range(0, max_abs_l + 1):
        total_l_tran = l * l_tran

        max_real_range_sqr = (length ** 2 - total_l_tran.imag() ** 2).upper()
        
        if max_real_range_sqr >= 0:
            max_real_range = RIF(max_real_range_sqr).sqrt()

            if l == 0:
                min_m = 1
            else:
                min_m = _min_int_in_interval(
                    (- total_l_tran.real() - max_real_range) / m_tran.real())
                
            max_m = _max_int_in_interval(
                    (- total_l_tran.real() + max_real_range) / m_tran.real())

            for m in range(min_m, max_m + 1):
                if gcd(m, l) == 1:
                    result.append((m,l))

    return result
                
def _max_int_in_interval(i):
    return i.upper().floor()

def _min_int_in_interval(i):
    return i.lower().ceil()


    
