"""

Giac [1] is included in Sage 6.8 (July 2015) and newer with Sage
providing a pexpect interpeter interface.  Additionally, there is a
Cython-based wrapper for Giac [2] with a Sage-specific incarnation [3]
which is installable as an optional Sage spkg::

  sage -i giacpy_sage

This module is to make it easy to access a version of Giac from within
Sage, preferring the Cython-based wrapper if available.

[1] https://www-fourier.ujf-grenoble.fr/~parisse/giac.html
[2] https://pypi.org/project/giacpy
[3] https://gitlab.math.univ-paris-diderot.fr/han/giacpy-sage

"""

from snappy.sage_helper import _within_sage

giac = None
have_giac = False
have_giacpy = False

if _within_sage:
    try:
        from giacpy_sage import libgiac as giac
        have_giac = have_giacpy = True
    except ImportError:
        # Older versions of giacpy_sage were called giacpy
        try:
            from giacpy import libgiac as giac
            have_giac = have_giacpy = True
        except ImportError:
            try:
                # Fallback to pexpect interface.
                from sage.all import giac
                have_giac = True
                have_giacpy = False
            except ImportError as e:
                raise ImportError(e.args[0] +
                                  'Probably your version of SageMath is too old '
                                  'to ship with Giac.')
