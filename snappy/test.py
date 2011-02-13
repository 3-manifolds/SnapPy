import doctest, inspect, os, sys, getopt
import snappy

optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
verbose = len(optlist) > 0
doctest.testmod(snappy.SnapPy, verbose=verbose)
