''' Twister is a program by Mark Bell, Tracy Hall and Saul Schleimer for 
constructing triangulations of surface bundles and Heegaard splittings from 
a description of a mapping class of a surface.

Twister is available from:
	https://github.com/MarkCBell/twister

See the included users guide ./docs/Twister.pdf for instructions for installing, 
testing and using Twister.

Provides:
  Surface - A class representing a squared surface with curves drawn on it.
  DT_drilling_surface - A function returning a standard surface from a knots 
    Dowker--Thistlethwaite code.
  DT_handles_surface - A function returning a standard surface from a knots 
    Dowker--Thistlethwaite code.
  surface_database - A set of all surfaces in Twister's surface database.
  version - A string containing the version number of the Twister kernel. '''

from .main import Surface, DT_drilling_surface, DT_handles_surface, surface_database, version
