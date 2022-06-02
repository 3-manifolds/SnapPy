import sys
import os

# python GLSLManifoldInsideView.py m004 "[0, 1, 0, 1, 0, 0, -1, -1]"

darwinTkMsg = """
On some versions of Mac OS X and Tk, it might be necessary to run the
following command to make the WASD navigation keys work properly:

    defaults write -g ApplePressAndHoldEnabled -bool false

The effect (disabling the ability to enter accented characters by, e.g.,
holding the e key) can be undone with:

    defaults write -g ApplePressAndHoldEnabled -bool true
"""

"""
An attempt at fixing the navigation keys (using pyobjc installed with pip):

    from Foundation import NSUserDefaults

    NSUserDefaults.standardUserDefaults().setBool_forKey_(False, 'ApplePressAndHoldEnabled')
    print(NSUserDefaults.standardUserDefaults().get('ApplePressAndHoldEnabled'))

"""

from snappy import Manifold
from snappy.gui import ViewerWindow

# Import raytracing directly from SnapPy source so that we can quickly
# iterate on shaders without the need to build/install SnapPy every time.

dev_path, dir_name = os.path.split(os.getcwd())
snappy_path, dir_name = os.path.split(dev_path)

sys.path.append(os.path.join(snappy_path, 'python'))

from raytracing.inside_viewer import *

def run_perf_test(): 
    gui = ViewerWindow(
        InsideViewer, Manifold("m004"))

    PerfTest(gui.view.widget)

def main(manifold, weights):
    if sys.platform == 'darwin':
        print(darwinTkMsg)

    gui = ViewerWindow(
        InsideViewer, manifold, weights = weights)
    gui.mainloop()

def to_index(s):
    var, face_num, tet_num = s.split('_')
    if var != 's':
        raise Exception("Not s")
    return 4 * int(tet_num) + int(face_num)

def check_weights(trig, weights):
    
    face_classes = trig._ptolemy_equations_identified_face_classes()
    pos_c2 = [ weights[to_index(face_class[2])] for face_class in face_classes ]
    neg_c2 = [ weights[to_index(face_class[3])] for face_class in face_classes ]

    for p, n in zip(pos_c2, neg_c2):
        if p != -n:
            raise Exception("Not matching")

    for row in trig._ptolemy_equations_boundary_map_2()[0]:
        if len(row) != len(pos_c2):
            raise Exception("Not matching")
        s = sum([e * p for e, p in zip(row, pos_c2)])
        if s != 0:
            raise Exception("Not in kernel")
    
if __name__ == '__main__':
    print(sys.argv)

    if sys.argv[1] == 'perf':
        run_perf_test()
    else:
        trig = Manifold(sys.argv[1])

        weights = None
        if len(sys.argv) == 3:
            weights = eval(sys.argv[2])

            check_weights(trig, weights)

        main(trig, weights)
