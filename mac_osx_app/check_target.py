"""
Tools for finding out which sdk was used to build an executable binary and
what its minimum deployment target is.
"""
from __future__ import print_function
from subprocess import Popen, PIPE
import sys

class TargetChecker:
    output_string = """
    {self.path} was built using the {self.sdk} sdk targeting {self.target}
    """
    
    def __init__(self, path):
        self.path = path
        self.target, self.sdk = self.get_target_and_sdk(path)

    def __repr__(self):
        return self.output_string.format(self=self)

    def get_target_and_sdk(self, path):
        proc = Popen(['otool', '-l', path], stdout=PIPE, stderr=PIPE)
        output, errors = proc.communicate()
        lines = iter(output.decode('ascii').split('\n'))
        for line in lines:
            parts = line.split()
            if len(parts) < 2:
                continue
            if parts[1] == 'LC_VERSION_MIN_MACOSX':
                break
        for line in lines:
            parts = line.split()
            if len(parts) < 2:
                continue
            if parts[0] == 'version':
                target = parts[1]
            elif parts[0] == 'sdk':
                sdk = parts[1]
                return target, sdk

class TkChecker(TargetChecker):
    Tk_path = '/Library/Frameworks/Tk.framework/Tk'
    Tcl_path = '/Library/Frameworks/Tcl.framework/Tcl'
    output_string = """
    Tk was built using the {self.Tk_sdk} sdk targeting {self.Tk_target}
    Tcl was built using the {self.Tcl_sdk} sdk targeting {self.Tcl_target}
    """
    
    def __init__(self):
        self.Tk_target, self.Tk_sdk = self.get_target_and_sdk(self.Tk_path)
        self.Tcl_target, self.Tcl_sdk = self.get_target_and_sdk(self.Tcl_path)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print(TkChecker())
    else:
        print(TargetChecker(sys.argv[1]))
