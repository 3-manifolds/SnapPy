"""

A self-signing notabot.cfg file looks like the following, where the
password is really all x's.

[developer]
username = person@somewhere.else
password = xxxx-xxxx-xxxx-xxxx
identity = -

[entitlements]
plist_file = entitlements.plist

[app]
app_path = dist/SnapPy.app
dmg_path = disk_images/SnapPy.dmg
"""

import os
import re
from math import ceil
from subprocess import check_call, call, Popen, PIPE
from notabot import Notarizer

def package_app(dmg_path):
    """
    Create a disk image containing the app, with a nice background and
    a symlink to the Applications folder.
    """
    image_dir, dmg_name = os.path.split(dmg_path)
    if not os.path.exists(image_dir):
        os.mkdir(image_dir)
    mount_name = "/Volumes/SnapPy"
    temp_path = dmg_path.replace('.dmg', '-tmp.dmg')
    # Make sure the dmg isn't currently mounted, or this won't work.  
    while os.path.exists(mount_name):
        print("Trying to eject " + mount_name)
        os.system("hdiutil detach " + mount_name)
    # Remove old dmgs if they exist.
    if os.path.exists(dmg_path):
        os.remove(dmg_path)
    if os.path.exists(temp_path):
        os.remove(temp_path)
    # Add symlink to /Applications if not there.
    if not os.path.exists("dist/Applications"):
        os.symlink("/Applications/", "dist/Applications")

    # Copy over the background and .DS_Store file.
    check_call(["rm", "-rf", "dist/.background"])
    os.mkdir("dist/.background")
    check_call(["cp", "dmg-maker/background.png", "dist/.background"])
    check_call(["cp", "dmg-maker/dotDS_Store", "dist/.DS_Store"])
        
    # Figure out the needed size.
    raw_size, errors = Popen(["du", "-sh", "dist"], stdout=PIPE).communicate()
    size, units = re.search("([0-9.]+)([KMG])", str(raw_size)).groups()
    new_size = "%d" % ceil(1.2 * float(size)) + units
    # Run hdiutil to build the dmg file.:
    check_call(["hdiutil", "makehybrid", "-hfs", "-hfs-volume-name", "SnapPy",
                "-hfs-openfolder", "dist", "dist", "-o", temp_path])
    check_call(["hdiutil", "convert", "-format", "UDZO", temp_path, "-o", dmg_path])
    os.remove(temp_path)
    # Delete the symlink to /Applications or egg_info will be glacial on newer setuptools.
    os.remove("dist/Applications")

class SnapPyNotarizer(Notarizer):

    def build_dmg(self):
        app_bundle = self.config['paths']['bundle_path']
        dmg_file = self.config['paths']['dmg_path']
        print('building dmg %s for %s'%(dmg_file, app_bundle))
        package_app(dmg_file)

if __name__ == '__main__':
    notarizer = SnapPyNotarizer('notabot.cfg')
    notarizer.run()
        
