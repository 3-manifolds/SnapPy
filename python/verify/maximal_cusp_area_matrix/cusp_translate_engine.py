from ...sage_helper import _within_sage

from ..upper_halfspace.finite_point import *

if _within_sage:
    from sage.all import vector, matrix

__all__ = ['CuspTranslateEngine']

class CuspTranslateEngine(object):
    def __init__(self, t0, t1):
        self.t0 = t0
        self.t1 = t1

        a = t0.real()
        b = t1.real()
        c = t0.imag()
        d = t1.imag()
        
        det = a * d - b * c

        self._matrix = matrix([[ d / det, -b / det],
                               [-c / det,  a / det]])

    def _to_vec(self, z):
        v = self._matrix * vector([z.real(), z.imag()])
        for e in v:
            if not (e.absolute_diameter() < 0.5):
                raise Exception("Too large interval")
        return v

    def _canonical_translates(self, z):
        def round_to_nearest_integers(interval):
            r = interval.round()
            return list(range(int(r.lower()), int(r.upper()) + 1))

        v = self._to_vec(z)
        integer_ranges = [ round_to_nearest_integers(i) for i in v ]
        
        for i in integer_ranges[0]:
            for j in integer_ranges[1]:
                yield z - i * self.t0 - j * self.t1

    def canonical_translates(self, finitePoint):
        """
        Test::

        sage: from sage.all import *
        sage: t0 = CIF(RIF(2.3, 2.30000000001), 3.4)
        sage: t1 = CIF(4.32, RIF(5.43, 5.4300000001))
        sage: c = CuspTranslateEngine(t0, t1)
        sage: z = CIF(0.23, 0.43)
        sage: t = RIF(5)
        sage: for i in range(-2, 3): # doctest: +NUMERIC6
        ...     for j in range(-2, 3):
        ...         print(list(c.canonical_translates(FinitePoint(z + i * t0 + j * t1, t))))
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.23000000000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.230000000000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.23000000000000001? + 0.43000000000000000?*I, 5)]
        [FinitePoint(0.230000000000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.23000000000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]
        [FinitePoint(0.2300000000? + 0.430000000?*I, 5)]

        sage: list(c.canonical_translates(FinitePoint(t0 / 2, t)))
        [FinitePoint(1.15000000000? + 1.7000000000000000?*I, 5), FinitePoint(-1.15000000000? - 1.7000000000000000?*I, 5)]

        sage: list(c.canonical_translates(FinitePoint(t1 / 2, t)))
        [FinitePoint(2.1600000000000002? + 2.7150000000?*I, 5), FinitePoint(-2.1600000000000002? - 2.7150000000?*I, 5)]

        sage: list(c.canonical_translates(FinitePoint(t0 / 2 + t1 / 2, t)))
        [FinitePoint(3.31000000001? + 4.4150000000?*I, 5), FinitePoint(-1.01000000000? - 1.015000000?*I, 5), FinitePoint(1.01000000000? + 1.0150000000?*I, 5), FinitePoint(-3.3100000000? - 4.415000000?*I, 5)]

        """

        for z in self._canonical_translates(finitePoint.z):
            yield FinitePoint(z, finitePoint.t)

    def _translate_to_match(self, z, targetZ):
        # to_vec checks that intervals can contain at most one integer.
        v = self._to_vec(targetZ - z)

        integers = [ interval.is_int()[1] for interval in v ]
        if None in integers:
            return None

        return z + integers[0] * self.t0 + integers[1] * self.t1

    def translate_to_match(self, finitePoint, targetFinitePoint):
        """

        sage: from sage.all import *
        sage: t0 = CIF(RIF(2.3, 2.30000000001), 3.4)
        sage: t1 = CIF(4.32, RIF(5.43, 5.4300000001))
        sage: c = CuspTranslateEngine(t0, t1)
        sage: z = CIF(RIF(0.23, 0.26), 0.43)
        sage: perturb = CIF(0.01, 0)
        sage: t = RIF(5)

        sage: for i in range(-2, 3): # doctest: +NUMERIC6
        ...     for j in range(-2, 3):
        ...         print(c.translate_to_match(FinitePoint(z + i * t0 + j * t1 + perturb, t), FinitePoint(z, t)))
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.43000000000000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        FinitePoint(0.3? + 0.430000000?*I, 5)
        
        sage: perturb = CIF(0.1, 0)
        sage: c.translate_to_match(FinitePoint(z + 2 * t0 + 1 * t1 + perturb, t), FinitePoint(z, t)) is None
        True
        
        
        """

        z = self._translate_to_match(finitePoint.z, targetFinitePoint.z)
        if z is None:
            return None
        return FinitePoint(z, finitePoint.t)
