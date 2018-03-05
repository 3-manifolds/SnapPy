#! /usr/bin/env python
"""
Creates a nice disk image, with background and /Applications symlink
for the app.

One issue here is that Snow Leopard uses a different (undocumented, of
course) format for the .DS_Store files than earlier versions, which makes
disk images created on it not work correctly on those systems.   Thus this "solution" uses a .DS_Store file created on Leopard as follows:

(1) Use Disk Utility to create a r/w DMG large enough to store everything and open it.

(2) Copy over the application and add a symlink to /Applications.

(3) Create a subdirectory ".background" containing the file "background.png".

(4) Open the disk image in the Finder and do View->Hide Tool Bar and then View->Show View Options.  To add the background picture inside the hidden directory, use cmd-shift-g in the file dialog.  Adjust everything to suit, close window and open it.   Then copy the .DS_Store file to dotDS_store.  

"""
import os, sys, re
from math import ceil

if sys.version_info.major == 2:
    name = "SnapPy"
else:
    name = "SnapPy3"
dist_dir = "../dist"
print('dmg name is %s'%name)

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
    os.system("hdiutil convert -format UDZO SnapPy-tmp.dmg -o %s.dmg"%name)
    os.remove("SnapPy-tmp.dmg")
    # Delete symlink to /Applications or egg_info will be glacial on newer setuptools.
    os.remove(dist_dir + "/Applications")
              
    
    
if __name__ == "__main__":
    main()



