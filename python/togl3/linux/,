# -*- tcl -*-
# Tcl package index file, version 1.1
#
if {[package vsatisfies [package provide Tcl] 9.0-]} {
    package ifneeded Togl 3.0 \
	    [list load [file join $dir libtcl9Togl3.0.so] [string totitle Togl]]
} else {
    package ifneeded Togl 3.0 \
	    [list load [file join $dir libTogl3.0.so] [string totitle Togl]]
}
