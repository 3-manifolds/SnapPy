"""
There are two Python interfaces to PHCpack, but both require
compiling from the Ada source, which is often a pain.  

This file implements a simple wrapper for the command-line version of
PHCpack.

>>> R = PolynomialRing(QQ, ['x', 'y'])
>>> I = R.ideal([R('x^2 + y^2 - 1'), R('(x -  1/2)^2 + y^2 - 1')])
>>> I.dimension()
0
>>> find_solutions(I, 2)[0]['coors']['y']
-0.9682458365518542212948163499456
"""

import re, os, tempfile
from sage.all import RealField, ComplexField, QQ, PolynomialRing


def ideal_to_file(ideal, filename):
    outfile = open(filename, 'w')
    polys = ideal.gens()    
    outfile.write('%d\n' % len(polys))
    for p in polys:
        outfile.write('   ' + repr(p) + ';\n')
    outfile.close()


def parse_file(filename, prec=53):
    RR = RealField(prec)
    CC = ComplexField(prec)
    data = open(filename).read()
    data = data.split('THE SOLUTIONS')[-1]
    ans = []
    solutions = re.findall('(solution \d+ : \s* start residual .*?) ==', data, re.DOTALL)
    for sol in solutions:
        kind = sol.split('=')[-1].strip()
        if kind == 'no solution':
            continue
        mult = int(re.search('^m : (\d+)', sol, re.MULTILINE).group(1))
        err = float(re.search('== err :\s+(.*?)\s+= ', sol).group(1))
        coors = re.findall('^ (.*) :\s+(\S*)\s+(\S*)', sol, re.MULTILINE)
        if kind.startswith('real'):
            coors = {var:RR(real) for var, real, imag in coors}
            ans.append({'kind':'real', 'mult':mult, 'err':err, 'coors':coors})
        elif kind.startswith('complex'):
            coors = {var:CC(RR(real), RR(imag)) for var, real, imag in coors}
            ans.append({'kind':'complex', 'mult':mult, 'err':err, 'coors':coors})

    num_real = int(re.search('Number of real solutions\s+:\s(.*).', data).group(1))
    num_complex = int(re.search('Number of complex solutions\s+:\s(.*).', data).group(1))
    kinds = [sol['kind'] for sol in ans]
    assert kinds.count('real') == num_real
    assert kinds.count('complex') == num_complex
    return ans

def find_solutions(ideal, doubles=1):
    assert doubles in [1, 2, 4]
    prec = 53*doubles
    tmpdir = tempfile.mkdtemp()
    infile = tmpdir + os.sep + 'polys.txt'
    outfile = tmpdir + os.sep + 'out.txt'    
    ideal_to_file(ideal, tmpdir + os.sep + 'polys.txt')
    flag = {1:'-b', 2:'-b2', 4:'-b4'}[doubles]
    os.system('phc ' + flag + ' ' + infile + ' ' + outfile)
    ans = parse_file(outfile, prec)
    os.system('rm -rf tmpdir')
    return ans

if __name__ == '__main__':
    import doctest
    doctest.testmod()
        

    
