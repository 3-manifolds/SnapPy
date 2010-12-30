#! /usr/bin/env python
"""
Creates a nice disk image, with background and /Applications symlink
for the app.

One issue here is that Snow Leopard uses a different (undocumented, of
course) format for the .DS_Store files than earlier versions, which makes
disk images created on it not work correctly on those systems.   Thus this
"solution" uses a .DS_Store file created on Leopard, and while it's a hack
at least it doesn't require the "AdiumApplescriptRunner" thing which gets
mess up (not its fault) by various Adobe products...
"""
import os, sys, re
from math import ceil

name = "SnapPy"
dist_dir = "../dist"

def main():
    # Make sure the dmg isn't currently mounted, or this won't work.  
    mount_name = "/Volumes/" + name
    while os.path.exists(mount_name):
        print("Trying to eject " + mount_name)
        os.system("hdiutil detach " + mount_name)
    # Remove old dmg if there is one
    while os.path.exists(name + ".dmg"):
        os.remove(name + ".dmg")
    while os.path.exists(name + "-tmp.dmg"):
        os.remove(name + "-tmp.dmg")
    # Add symlink to /Applications if not there:
    if not os.path.exists(dist_dir + "/Applications"):
        os.symlink("/Applications/", dist_dir + "/Applications")

    # copy over the background and .DS_Store file
    os.system("rm -Rf " + dist_dir + "/.background")
    os.system("mkdir " + dist_dir + "/.background")
    os.system("cp background.png " + dist_dir + "/.background")
    os.system("cp dotDS_Store " + dist_dir + "/.DS_Store")
    
    # figure out the needed size:
    raw_size = os.popen("du -sh " + dist_dir).read()
    size, units = re.search("([0-9.]+)([KMG])", raw_size).groups()
    new_size = "%d" % ceil(1.2 * float(size)) + units
    # Run the main script:
    os.system("hdiutil makehybrid -hfs -hfs-volume-name SnapPy -hfs-openfolder %s %s -o SnapPy-tmp.dmg" % (dist_dir, dist_dir))
    os.system("hdiutil convert -format UDZO SnapPy-tmp.dmg -o SnapPy.dmg")
    os.remove("SnapPy-tmp.dmg")
              
    
    
if __name__ == "__main__":
    main()



