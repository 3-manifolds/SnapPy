#
# Tcl package index file
#
package ifneeded mactoolbar 1.0 "
    package require Tk 8.5-
    if {\"AppKit\" ni \[winfo server .\]} {error {TkAqua Cocoa required}}
        load [list ][file join $dir libmactoolbar1.0.dylib] mactoolbar"
