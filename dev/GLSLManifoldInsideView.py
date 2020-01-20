from __future__ import print_function

"""
Cheats:

defaults write -g ApplePressAndHoldEnabled -bool false



"""

from raytracing.manifold_inside_view import *
from snappy import Manifold

def run_perf_test(): 
    gui = InsideManifoldGUI(Manifold("m004"))

    PerfTest(gui.main_widget)

def main(manifold):
    gui = InsideManifoldGUI(manifold)
    gui.main_widget.focus_set()
    gui.window.mainloop()
    
if __name__ == '__main__':
    print(sys.argv)

    if sys.argv[1] == 'perf':
        run_perf_test()
    else:
        main(Manifold(sys.argv[1]))
