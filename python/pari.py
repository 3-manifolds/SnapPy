"""
Import pari and associated classes and functions here, to be more DRY,
while supporting both old and new versions of cypari and sage.pari and
accounting for all of the various idiosyncrasies.
"""
from pkg_resources import parse_version
from .sage_helper import _within_sage

if _within_sage:
    from sage.version import version as sage_version
    sage_version = parse_version(sage_version)
    if sage_version < parse_version('7.5'):
        from sage.libs.pari.gen import gen as Gen
        if sage_version < parse_version('6.1'):
            from sage.libs.pari.gen import pari
            from sage.libs.pari.gen import (prec_words_to_dec,
                                            prec_words_to_bits,
                                            prec_bits_to_dec,
                                            prec_dec_to_bits)
        else:
            from sage.libs.pari.pari_instance import pari
            from sage.libs.pari.pari_instance import (prec_words_to_dec,
                                                      prec_words_to_bits,
                                                      prec_bits_to_dec,
                                                      prec_dec_to_bits)
    elif sage_version < parse_version('7.6'):
        from sage.libs.cypari2 import pari
        from sage.libs.cypari2.gen import gen as Gen
        from sage.libs.cypari2.pari_instance import (
            prec_words_to_dec,
            prec_words_to_bits,
            prec_bits_to_dec,
            prec_dec_to_bits)
    elif sage_version < parse_version('8.0'):
        from sage.libs.pari import pari
        from sage.libs.cypari2 import Gen
        from sage.libs.cypari2.pari_instance import (
            prec_words_to_dec,
            prec_words_to_bits,
            prec_bits_to_dec,
            prec_dec_to_bits)
    else: # Sage 8.0 or newer
        from sage.libs.pari import pari
        from cypari2 import Gen
        from cypari2.pari_instance import (
            prec_words_to_dec,
            prec_words_to_bits,
            prec_bits_to_dec,
            prec_dec_to_bits)
        
    from sage.all import PariError
    shut_up  = lambda : None
    speak_up = lambda : None   
    
else: # Plain Python, use CyPari
    import cypari
    cypari_version = parse_version(cypari.__version__)
    try:
        from cypari import pari
    except ImportError:
        from cypari.gen import pari
    if cypari_version < parse_version('2.2.0'):
        from cypari.gen import (
            gen as Gen,
            PariError,
            prec_words_to_dec,
            prec_words_to_bits,
            prec_bits_to_dec,
            prec_dec_to_bits)
    else:
        from cypari._pari import (
            Gen,
            PariError,
            prec_words_to_dec,
            prec_words_to_bits,
            prec_bits_to_dec,
            prec_dec_to_bits)
    shut_up  = lambda : pari.shut_up()
    speak_up = lambda : pari.speak_up()

