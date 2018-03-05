# Has been fixed!

import snappy
sample_isosig = 'hLMzMkbcdefggghhhhhhfo'   # v0000
sample_file_contents = snappy.Manifold(sample_isosig)._to_string()

def show_progress(i):
    if i % 1000 == 0:
        print(i)
    return i + 1


def from_file():
    """
    No leak here!
    """
    i = 0
    while True:
        M = snappy.Manifold(sample_file_contents)
        i = show_progress(i)
        del M
        snappy.SnapPy.check_SnapPea_memory()

def from_isosig():
    """
    Leaks about 1Gb/minute on NMD's OS X laptop
    """
    i = 0
    while True:
        M = snappy.Manifold(sample_isosig)
        i = show_progress(i)
        del M
        snappy.SnapPy.check_SnapPea_memory()

def to_isosig():
    """
    Leaks about 0.4Gb/minute on NMD's OS X laptop
    """
    i = 0
    while True:
        M = snappy.Manifold(sample_file_contents)
        M.triangulation_isosig()
        i = show_progress(i)
        del M
        snappy.SnapPy.check_SnapPea_memory()

from_isosig()
