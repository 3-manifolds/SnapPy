"branch","m004","m004 - 10","m003(-3,1)","m003(-3,1) - 10","o9_00637","m129(0,0)(-55,34)","m129(0,0)(-55,34) - 7","gLLMQcbeefefpnamups(2,1)","gLLMQcbeefefpnamups(2,1) - 7"

from snappy import Manifold
import time
import sys


class Timer:
    def __init__(self, times):
        self.times = times
    
    def __enter__(self):
        self.start_time = time.time()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end_time = time.time()
        self.elapsed_time = self.end_time - self.start_time
        self.times.append(self.elapsed_time)

def assert_real_length(geodesic, real_length, epsilon=1e-9):
    geodesic_real_length = geodesic.length.real()
    RF = geodesic_real_length.parent()
    if not abs(geodesic_real_length - RF(real_length)) < RF(epsilon):
        raise AssertionError(
            "Expected length %r. Got: %r." % (
                real_length, geodesic_real_length))

def assert_real_lengths(geodesics, real_lengths, epsilon=1e-9):
    if len(geodesics) != len(real_lengths):
        raise AssertionError(
            "Expected %d geodesics. Got %d." % (
                len(real_lengths), len(geodesics)))
    for geodesic, real_length in zip(geodesics, real_lengths):
        assert_real_length(geodesic, real_length, epsilon)

branch = sys.argv[1]

times = []

verified = True
bits_prec = 100
###############################################################################
M = Manifold('m004')

############################################################
with Timer(times):
    systole = next(M.length_spectrum_alt_gen(verified=verified, bits_prec=bits_prec))
assert_real_length(systole, 1.08707014499573909978527280123)

############################################################
with Timer(times):
    geodesics = M.length_spectrum_alt(count=10, verified=verified, bits_prec=bits_prec)
assert_real_lengths(
    geodesics,
    [1.0870701449957390997852728012,
     1.0870701449957390997852728012,
     1.6628858910586210756524850391,
     1.6628858910586210756524850391,
     1.6628858910586210756524850391,
     1.6628858910586210756524850391,
     1.7251092553241220643713855890,
     1.7251092553241220643713855889,
     2.1741402899914781995705456025,
     2.1741402899914781995705456025,
     2.1741402899914781995705456025,
     2.1741402899914781995705456025])

###############################################################################
M = Manifold('m003(-3,1)')

############################################################
with Timer(times):
    systole = next(M.length_spectrum_alt_gen(verified=verified, bits_prec=bits_prec))
assert_real_length(systole, 0.584603685017986969321)

bits_prec=1000

############################################################
with Timer(times):
    geodesics = M.length_spectrum_alt(count=10, verified=verified, bits_prec=bits_prec)
assert_real_lengths(
    geodesics,
    [0.584603685017987,
     0.584603685017987,
     0.584603685017987,
     0.7941346629918657,
     0.7941346629918657,
     0.7941346629918657,
     1.2898511615162902,
     1.2898511615162902,
     1.2898511615162902,
     2.1791501442612917,
     2.1791501442612917,
     2.1791501442612917])

###############################################################################
M = Manifold('o9_00637')
bits_prec = 200
verified = False

############################################################
with Timer(times):
    systole = next(M.length_spectrum_alt_gen(verified=verified, bits_prec=bits_prec))
assert_real_length(systole, 0.00164352273449)

###############################################################################
M = Manifold('m129(0,0)(-55,34)')
bits_prec = 200
verified = True

############################################################
with Timer(times):
    systole = next(M.length_spectrum_alt_gen(verified=verified, bits_prec=bits_prec))
assert_real_length(systole, 0.00164352273449)

############################################################
with Timer(times):
    geodesics = M.length_spectrum_alt(count=7, verified=verified, bits_prec=bits_prec)
assert_real_lengths(
    geodesics,
    [0.0016435227344936082,
     1.061898381294215,
     1.061898381294215,
     1.7625566130602213,
     1.7625566130602213,
     1.764532805884113,
     1.764532805884113])

###############################################################################
M = Manifold('gLLMQcbeefefpnamups(2,1)')

############################################################
with Timer(times):
    systole = next(M.length_spectrum_alt_gen(verified=verified, bits_prec=bits_prec))
assert_real_length(systole, 0.83144294552931)

############################################################
with Timer(times):
    geodesics = M.length_spectrum_alt(count=7, verified=verified, bits_prec=bits_prec)
assert_real_lengths(
    geodesics,
    [0.83144294552931,
     0.83144294552931,
     0.83144294552931,
     0.86255462766206,
     0.86255462766206,
     0.86255462766206,
     1.31695789692481,
     1.31695789692481,
     1.31695789692481])

print('"%s",%s' % (branch, ','.join('%.4f' % time for time in times)))

# v3210(-3,1)
# v2923(-1,4)
