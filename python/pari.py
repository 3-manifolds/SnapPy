"""
Import pari and associated classes and functions here, to be more DRY,
while supporting both old and new versions of cypari and sage.pari and
accounting for all of the various idiosyncrasies.
"""
try: # Sage
    from sage.libs.pari.gen import gen
    try:
        from sage.libs.pari.gen import pari
        from sage.libs.pari.gen import (prec_words_to_dec,
                                        prec_words_to_bits,
                                        prec_bits_to_dec,
                                        prec_dec_to_bits)
    except ImportError: # Sage 6.1 or later:
        from sage.libs.pari.pari_instance import pari
        from sage.libs.pari.pari_instance import (prec_words_to_dec,
                                                  prec_words_to_bits,
                                                  prec_bits_to_dec,
                                                  prec_dec_to_bits)
    from sage.all import PariError
    shut_up  = lambda : None
    speak_up = lambda : None   
    
except ImportError: # CyPari
    try:
        from cypari.gen import pari
        from cypari.gen import (
            gen, PariError,
            prec_words_to_dec,
            prec_words_to_bits,
            prec_bits_to_dec,
            prec_dec_to_bits)
        shut_up  = lambda : pari.shut_up()
        speak_up = lambda : pari.speak_up()
    except ImportError: # CyPari Version 2:
        from cypari import pari
        from cypari.gen import gen
        from cypari.pari_instance import (
            prec_words_to_dec,
            prec_words_to_bits,
            prec_bits_to_dec,
            prec_dec_to_bits)
        shut_up  = lambda : None
        speak_up = lambda : None
    from cypari.gen import PariError
