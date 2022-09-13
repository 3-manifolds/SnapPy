"""
Import pari and associated classes and functions here, to be more DRY,
while supporting both old and new versions of cypari and sage.pari and
accounting for all of the various idiosyncrasies.
"""
from pkg_resources import parse_version
from .sage_helper import _within_sage

if _within_sage:
    from sage.version import version as s_version
    sage_version = parse_version(s_version)
    if sage_version < parse_version('9.0'):
        raise ValueError("you need a more recent version of SageMath")
    from sage.libs.pari import pari
    from cypari2 import Gen
    from cypari2.pari_instance import (prec_words_to_dec,
                                       prec_words_to_bits,
                                       prec_bits_to_dec,
                                       prec_dec_to_bits)
    from sage.all import PariError
    shut_up = lambda: None
    speak_up = lambda: None

else:  # Plain Python, use CyPari
    import cypari
    cypari_version = parse_version(cypari.__version__)
    if cypari_version < parse_version('2.3'):
        raise ValueError("you need a more recent version of CyPari")
    from cypari import pari
    from cypari._pari import (Gen,
                              PariError,
                              prec_words_to_dec,
                              prec_words_to_bits,
                              prec_bits_to_dec,
                              prec_dec_to_bits)
    shut_up = lambda: pari.shut_up()
    speak_up = lambda: pari.speak_up()
