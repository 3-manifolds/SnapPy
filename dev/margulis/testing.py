from supremal_margulis_number import margulis_number
from predicate import is_margulis_number
import time
import sys
from snappy import OrientableCuspedCensus
from snappy import Manifold

def test_examples():
    M = Manifold("m125")
    print(margulis_number(M))
    M = Manifold("m004")
    print(margulis_number(M))
    M = Manifold("o9_10000")
    print(margulis_number(M))
    print(margulis_number(M, bits_prec=500, verified=True))

    M = Manifold('m129(5,4)') # Same as m076, involves core curve
    print(margulis_number(M))
    M = Manifold("m076")
    print(margulis_number(M))

def do_census():

    for i in range(int(sys.argv[1]) * 3000, min(61911, (1 + int(sys.argv[1])) * 3000)):
        M = OrientableCuspedCensus[i]

        sys.stdout.flush()
        s = time.process_time()
        for bits_prec in [500, 1000, 2000]:
            try:
                epsilon, thin_part, collisions = margulis_number(
                    M, bits_prec=bits_prec, verified=True, include_thin_part=True)
                break
            except Exception as e:
                exception = e
        else:
            raise exception
        time_epsilon = time.process_time() - s

        if epsilon.is_NaN():
            raise Exception("NaN")

        if not epsilon > 0:
            raise Exception("Zero epsilon")

        if epsilon.diameter() > 1e-10:
            raise Exception("Big diameter")

        RF = epsilon.parent()
        delta = RF(1e-25)

        s = time.process_time()

        for bits_prec in [500, 1000, 2000]:
            try:
                m1 = is_margulis_number(M, epsilon - delta, bits_prec=bits_prec, verified=True)[0]
                break
            except Exception as e:
                e = e
        if not m1:
            raise Exception("Margulis 1")

        for bits_prec in [500, 1000, 2000]:
            try:
                m2 = is_margulis_number(M, epsilon + delta, bits_prec=bits_prec, verified=True)[0]
                break
            except Exception as e:
                e = e
        if m2:
            raise Exception("Margulis 1")

        time_predicate = time.process_time() - s

        print('%s,%r,%f,%f,"%r","%r"' % (
            M.name(), epsilon, time_epsilon, time_predicate, thin_part, collisions))
        sys.stdout.flush()

def do_closed_census():
    for i in range(int(sys.argv[1]) * 2000, min(11031, (1 + int(sys.argv[1])) * 2000)):
        # Fails on Closed114.tri

        name = 'Closed%d.tri' % (i + 1)
        M = Manifold("ClosedManifolds/%s" % name)

        if False and name == 'Closed3.tri':
            # Three-fold cover is
            # M=Manifold('gLLMQcbeefefpjaqupw_BBab(1,0)')
            # gLLMQcbeefefpjaqupw(1,1)
            # m247(-1,3)

            print("%s" % name)
            continue

        sys.stdout.flush()
        s = time.process_time()
        for bits_prec in [500, 1000, 2000, 5000]:
            try:
                epsilon, thin_part, collisions = margulis_number(
                    M, bits_prec=bits_prec, verified=True, include_thin_part=True)
                break
            except Exception as e:
                exception = e
        else:
            print('%s,,,,"%s"' % (name, str(exception).replace('\n',' ')))
            continue
        time_epsilon = time.process_time() - s

        if epsilon.is_NaN():
            raise Exception("NaN")

        if not epsilon > 0:
            raise Exception("Zero epsilon")

        if epsilon.diameter() > 1e-10:
            raise Exception("Big diameter")

        print('%s,%r,%f,"%r","%r"' % (
            name, epsilon, time_epsilon, thin_part, collisions))
        sys.stdout.flush()


if __name__ == '__main__':
    do_census()
    #do_closed_census()
    # test_examples()
