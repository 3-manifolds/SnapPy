import sys
import IPython
from snappy import *
from IPython.Shell import IPShellEmbed
from snappy.phone_home import needs_updating

try:
    from site import _Printer

    copyright = _Printer(name='copyright', data = str(copyright) + """

SnapPy: Copyright (c) 2009-present by Marc Culler, Nathan Dunfield and others.
        All Rights Reserved
""")

    credits = _Printer(name='credits', data = """
Python: 
    Thanks to CWI, CNRI, BeOpen.com, Zope Corporation and a cast of thousands
    for supporting Python development.  See www.python.org for more information.

IPython:
    Fernando Perez, Janko Hauser, Nathan Gray, and many users.
    See http://ipython.scipy.org for more information.

SnapPy:
    Marc Culler, Nathan Dunfield, Jeff Weeks, and many topologists.
    See http://snappy.computop.org for more information.
""")

    license = _Printer(name='license', data = """
Type license() to see the full text of the Python license.
SnapPy is distributed under the GNU Public License.
See http://www.gnu.org/copyleft/gpl.html
""")

except:
    pass


SnapPy_banner = """
    Hi.  It's SnapPy.
    SnapPy is based on the SnapPea kernel, written by Jeff Weeks.
    Type Manifold? to get started.
    Type "copyright", "credits", or "license" for more information."""

status = needs_updating()
if status:
    SnapPy_banner += "\n    **Please upgrade to %s from %s via http://snappy.computop.org**" % status

def SnapPy_showtraceback(exc_tuple = None,filename=None,tb_offset=None):
    if exc_tuple is None and the_shell.IP.tracebacks == False:
        type, value, tb = sys.exc_info()
        print '\033[0;31m%s:\033[0m  %s'%(type.__name__, value)
    else:
        the_shell.IP.IP_showtraceback(exc_tuple, filename, tb_offset)

# Instantiate an IPython shell
the_shell = IPShellEmbed(banner=SnapPy_banner)

# Tinker with its tracebacks
the_shell.IP.tracebacks = False
the_shell.IP.IP_showtraceback = the_shell.IP.showtraceback
the_shell.IP.showtraceback = SnapPy_showtraceback
# To restore tracebacks: __IP.tracebacks = True

# Start it up
if __name__ == "__main__":
    the_shell()
