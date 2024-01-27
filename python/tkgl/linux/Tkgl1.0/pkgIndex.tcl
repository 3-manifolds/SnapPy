# -*- tcl -*-
# Tcl package index file, version 1.1
#
if {[package vsatisfies [package provide Tcl] 9.0-]} {
    package ifneeded Tkgl 1.0 \
	    [list load [file join $dir libtcl9Tkgl1.0.so] [string totitle Tkgl]]
} else {
    package ifneeded Tkgl 1.0 \
	    [list load [file join $dir libTkgl1.0.so] [string totitle Tkgl]]
}
