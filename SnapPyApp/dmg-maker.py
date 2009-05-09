#! /usr/bin/env python
"""

Utility for creating compressed MacOS disk images of the given
directory with the following special feature.

* When mounted, they pop up a window with the disk contents.

Usage:

    dmg-maker.py   directory   image_name.dmg
"""

import os, sys, re, time

def create_image(directory, image_name):
    """
    Creates a compressed disk image containing
    the contents of directory
    """

    # Remove old file, if any
    
    os.system("rm " + image_name)

    # Create new image
    
    os.system("hdiutil create " + 
              " -volname " + image_name[:-4] + 
              " -srcfolder " + directory + 
              " -format UDBZ "+
              " -o " + image_name
              )

def make_image_auto_open_window(image_name):
    #os.system("hdiutil attach -readwrite -noverify -noautoopen " + image_name)
    #os.system("sudo bless --folder " + vol_name)
    #os.system("sudo bless --openfolder " + vol_name)
    # os.system("umount " + vol_name)
    pass 

if __name__ == "__main__":
    directory, image_name = sys.argv[1:]
    assert len(sys.argv) == 3 and image_name[-4:] == ".dmg"
    create_image(directory, image_name)
    
    
