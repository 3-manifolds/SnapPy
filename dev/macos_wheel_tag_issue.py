import os, sys
from wheel.macosx_libfile import extract_macosx_min_system_version

for (dirpath, dirnames, filenames) in os.walk(sys.argv[1]):
    for filename in filenames:
        if filename.endswith('.dylib') or filename.endswith('.so'):
            lib_path = os.path.join(dirpath, filename)
            min_ver = extract_macosx_min_system_version(lib_path)
            if min_ver is not None:
                print(min_ver, lib_path)

