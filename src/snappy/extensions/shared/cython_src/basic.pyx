def valid_index(i, n, format_str):
    """
    Return range(n)[i] or raises a nicely formatted IndexError
    using format_str.

    This does several things for us::

        * avoid that cython happily converts a float to an int, so a call
          such as M.dehn_fill((1,0), 0.6) would succeed.
        * converts a negative index such as -1 to n - 1
        * checks that i is in the right range
        * supports Sage and numpy Integers: they are not python int's but
          have __index__ (see PEP 357)

    It is faster than reimplementing these behaviors.
    """
    try:
        return range(n)[i]
    except IndexError:
        raise IndexError(format_str % i)
