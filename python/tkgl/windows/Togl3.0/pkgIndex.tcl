if {[package vsatisfies [package provide Tcl] 9.0-]} { 
package ifneeded Togl 3.0 [list load [file join $dir tcl9Togl30t.dll]] 
} else { 
package ifneeded Togl 3.0 [list load [file join $dir Togl30t.dll]] 
} 
