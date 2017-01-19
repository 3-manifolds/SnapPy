#! /usr/bin/env python

docs = """
Add a "release candidate" tag to wheels/eggs/tarballs.  
"""
   

import sys, os, re, argparse
parser = argparse.ArgumentParser(description='adds rc tags for uploading to testpypi.',
                                 epilog=docs, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-r', '--rctag', help='index for rc tag to be appended to version (e.g. -r2 -> rc2)',
                    action='store', required=True)
parser.add_argument('files', nargs='+')

suffixes = re.compile('(.tar.gz|.zip|.whl|.egg)$')
version_tag = re.compile('([^-]+)-([^-]+?)(rc\d+|)(-|.tar.gz|.zip)')
    
if __name__ == '__main__':
    args = parser.parse_args()
    for file in args.files:
        if suffixes.search(file) is None:
            print('Skipping %s' % file)
        name = os.path.basename(file)
        path = os.path.dirname(file)
        tag = 'rc' + args.rctag 
        new_name = version_tag.sub('\g<1>-\g<2>%s\g<4>'%tag, name, 1)
        print('%s --> %s' % (name, new_name))
        os.rename(os.path.join(path, name), os.path.join(path, new_name))
