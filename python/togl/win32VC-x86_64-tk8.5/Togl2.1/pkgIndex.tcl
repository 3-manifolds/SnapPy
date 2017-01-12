#
# Tcl package index file
#
package ifneeded Togl 2.1 \
    [list load [file join $dir Togl21.dll]]
package ifneeded Togl 2.1 [list load [file join $dir Togl21.dll] Togl] 
