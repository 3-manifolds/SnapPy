"""
IMPORTANT: Python only recognises this as a doc string if there is
nothing before it. In particular, add any includes after the doc string.

Doc tests that we do not want to show up in the documentation.

isometry_signature
------------------

>>> M = Manifold('m004(1,2)')
>>> M.isometry_signature()
'cPcbbbiht(3,2)'
sage: M.isometry_signature(verified=True)
'cPcbbbiht(3,2)'

>>> M = ManifoldHP('m004(1,2)')
>>> M.isometry_signature()
'cPcbbbiht(3,2)'
sage: M.isometry_signature(verified=True) # Test bug reported by Nathan
'cPcbbbiht(3,2)'

# Cases where the drilled manifold had a canonical cell decomposition
# with non-tetrahedral cells.

>>> M = Manifold('m137(3,2)')
>>> M.isometry_signature()
'sLLvwzvQPAQPQccghmiljkpmqnoorqrrqfafaoaqoofaoooqqaf(3,2)'
sage: M.isometry_signature(verified=True)
'sLLvwzvQPAQPQccghmiljkpmqnoorqrrqfafaoaqoofaoooqqaf(3,2)'

>>> M = ManifoldHP('m137(3,2)')
>>> M.isometry_signature()
'sLLvwzvQPAQPQccghmiljkpmqnoorqrrqfafaoaqoofaoooqqaf(3,2)'
sage: M.isometry_signature(verified=True)
'sLLvwzvQPAQPQccghmiljkpmqnoorqrrqfafaoaqoofaoooqqaf(3,2)'

"""

if not __doc__:
    raise Exception("doc string with tests was not recognized.")

from . import Manifold, ManifoldHP
