import snappy
import time

iso16 = 'qLLvLLzAPQQkcehjlnhnopommpponiiaciijpxlanfggfk_baBaDBaBbB'

iso20 = 'uLLLLvPzMLAPQQccefemllnkmqorptstsrtsiitditpatgevpcmppuggn_abBaBbHg'

iso34 = 'ILLALLMvvAzLLPPMMvzAQQQkbcdeghiklqoqstwxzvAuxCGBGAEHCDHGFFHtsfxjxajmdobxadrfcdgoaamlloocnkkoks_BaCB'

iso48 = 'WLLLwvzvAAQALzPMQLzLvAPzQMzQzAPLQcdfemopklnqmnntsxwuywxzxACIDGHIHLHKMOMPPPQRTSVUVViceaiocvfvvwfnaamofdvpaikdcrvfbbhopcafaacaaaaaagb_baDB'

iso66 = '-ccbLLLvPAvvQQvvLPPzPLMMLMwvAzAzQvPwLQPLPPAzQPQQcaeafaiajaiakanalalamaqaqasauaxazaAayaEaBaFazaEaHaBaGaJaIaEaOaLaLaKaRaVaUaVaVaQaZaYa3a2a4aWaXaXa6a1a0a1a9a8a5a7a-a-a+abbab+a9a9aabbbbbiiegvgutnfksktaafffaaabaaabffgaqevadwqucaapaaaxdifocaavpiqqdjrcrgea_bBcB'

iso90 = '-cAbLLLLLLLvzvvvPAvLPvAvAAMvzwMwzvLLQzvQAQQQQzwwLPQvMPQQMPQQLQQQcadafahaiamaiapaqaDaravauawaxaKaCaRaBaBaBaMaMaWaIaYa2a0aZa2a9a4aPaQaQaSaVaUaVaabhbgbcb1a1aYa1aYaabbbbbab5a8a6akblbqbebfbdbibkbvbobsbrbsbrbrbtbtbubxbmbqbnbmbobsbpbqbzbybxbxbvbwbwbzbzbicegvivfiuatboftfmsgnqqovadaiqalbkbbcvcmapadqmthskgsvfglmajbfqimiqtqdplxmdvggvjrgsaaeefcwfo_bBRqp'

class DummyNumber:
    def __init__(self, val):
        self.val = val

    def __mul__(self, other):
        if isinstance(other, DummyNumber):
            return DummyNumber(self.val * other.val)

    def __rmul__(self, other):
        if isinstance(other, int):
            return DummyNumber(other * self.val)

    def __pow__(self, other):
        if isinstance(other, int):
            return DummyNumber(self.val ** other)

    def __sub__(self, other):
        if isinstance(other, int):
            return DummyNumber(self.val - other)

    def __rsub__(self, other):
        if isinstance(other, int):
            return DummyNumber(other - self.val)

    def __repr__(self):
        return repr(self.val)


class DummyNumber2:
    def __init__(self, val):
        self.val = val

    def _binop(self, operator, other):
        operand = other.val if isinstance(other, DummyNumber2) else other
        return DummyNumber2(operator(operand))

    def __mul__(self, other):
        return self._binop(self.val.__mul__, other)

    def __rmul__(self, other):
        return self._binop(self.val.__rmul__, other)

    def __pow__(self, other):
        return self._binop(self.val.__pow__, other)

    def __sub__(self, other):
        return self._binop(self.val.__sub__, other)

    def __rsub__(self, other):
        return self._binop(self.val.__rsub__, other)

    def __repr__(self):
        return repr(self.val)


try:
    from snappy.number import bit_precision
except ImportError:
    pass

class DummyNumber3:
    def __init__(self, val, prec=53):
        self.val = val
        self.prec = prec

    def _binop(self, operator, other):
        if isinstance(other, DummyNumber3):
            prec = min(self.prec, other.prec)
            with bit_precision(prec):
                return DummyNumber3(operator(other.val), prec)
        with bit_precision(self.prec):
            return DummyNumber3(operator(other), self.prec)

    def __mul__(self, other):
        return self._binop(self.val.__mul__, other)

    def __rmul__(self, other):
        return self._binop(self.val.__rmul__, other)

    def __pow__(self, other):
        return self._binop(self.val.__pow__, other)

    def __sub__(self, other):
        return self._binop(self.val.__sub__, other)

    def __rsub__(self, other):
        return self._binop(self.val.__rsub__, other)

    def __repr__(self):
        return repr(self.val)



def eval_gluing_equation(eqn, acb_shapes):
    """
    Evaluate the product of cross ratios in an edge equation.
    The result will be 1 if the equation is satisfied.
    """
    a, b, c = eqn
    ans = c
    for i, z in enumerate(acb_shapes):
        ans = ans * (z**int(a[i]) * (1 - z)**int(b[i]))
    return ans


def gluing_equation_errors(eqns, acb_shapes):
    """ A list containing the difference from 1 for each equation."""
    return [eval_gluing_equation(eqn, acb_shapes) - 1 for eqn in eqns]


def main_number():
    total_time = 0
    for iso in [iso16, iso20, iso34, iso48, iso66, iso90]:
        M = snappy.Manifold(iso)
        eqns = M.gluing_equations('rect')
        shapes = M.tetrahedra_shapes('rect')
        start = time.time()
        for i in range(100):
            gluing_equation_errors(eqns, shapes)
        total_time += (time.time() - start)

    print('number', total_time)


def main_complex():
    total_time = 0
    for iso in [iso16, iso20, iso34, iso48, iso66, iso90]:
        M = snappy.Manifold(iso)
        eqns = M.gluing_equations('rect')
        shapes = [complex(z) for z in M.tetrahedra_shapes('rect')]
        start = time.time()
        for i in range(100):
            gluing_equation_errors(eqns, shapes)
        total_time += (time.time() - start)

    print('complex', total_time)


def main_pari():
    total_time = 0
    for iso in [iso16, iso20, iso34, iso48, iso66, iso90]:
        M = snappy.Manifold(iso)
        eqns = M.gluing_equations('rect')
        shapes = [z.gen for z in M.tetrahedra_shapes('rect')]
        start = time.time()
        for i in range(100):
            gluing_equation_errors(eqns, shapes)
        total_time += (time.time() - start)

    print('bare pari', total_time)


def main_flint():
    total_time = 0
    for iso in [iso16, iso20, iso34, iso48, iso66, iso90]:
        M = snappy.Manifold(iso)
        eqns = M.gluing_equations('rect')
        shapes = [z.flint_obj for z in M.tetrahedra_shapes('rect')]
        start = time.time()
        for i in range(100):
            gluing_equation_errors(eqns, shapes)
        total_time += (time.time() - start)

    print('bare flint', total_time)


def main_flint_and_dummy():
    total_time = 0
    for iso in [iso16, iso20, iso34, iso48, iso66, iso90]:
        M = snappy.Manifold(iso)
        eqns = M.gluing_equations('rect')
        shapes = [DummyNumber(z.flint_obj) for z in M.tetrahedra_shapes('rect')]
        start = time.time()
        for i in range(100):
            gluing_equation_errors(eqns, shapes)
        total_time += (time.time() - start)

    print('dummy num + flint', total_time)


def main_flint_and_dummy2():
    total_time = 0
    for iso in [iso16, iso20, iso34, iso48, iso66, iso90]:
        M = snappy.Manifold(iso)
        eqns = M.gluing_equations('rect')
        shapes = [DummyNumber2(z.flint_obj) for z in M.tetrahedra_shapes('rect')]
        start = time.time()
        for i in range(100):
            gluing_equation_errors(eqns, shapes)
        total_time += (time.time() - start)

    print('dummy2 num + flint', total_time)


def main_flint_and_dummy3():
    total_time = 0
    for iso in [iso16, iso20, iso34, iso48, iso66, iso90]:
        M = snappy.Manifold(iso)
        eqns = M.gluing_equations('rect')
        shapes = [DummyNumber3(z.flint_obj) for z in M.tetrahedra_shapes('rect')]
        start = time.time()
        for i in range(100):
            gluing_equation_errors(eqns, shapes)
        total_time += (time.time() - start)

    print('dummy3 num + flint', total_time)


def main_complex_and_dummy():
    total_time = 0
    for iso in [iso16, iso20, iso34, iso48, iso66, iso90]:
        M = snappy.Manifold(iso)
        eqns = M.gluing_equations('rect')
        shapes = [DummyNumber(complex(z)) for z in M.tetrahedra_shapes('rect')]
        start = time.time()
        for i in range(100):
            gluing_equation_errors(eqns, shapes)
        total_time += (time.time() - start)

    print('dummy num + flint', total_time)

#main_number()
#main_complex()
#main_pari()
#main_flint()
#main_flint_and_dummy()
#main_flint_and_dummy2()
#main_flint_and_dummy3()
#main_complex_and_dummy()
