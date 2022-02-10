#e.g. python3 get_pol.py mpsolve_pol/src/tests/unisolve/chrma22.pol
import os
import sys
import mpmath as mp
import shlex
import subprocess
import re
import argparse

#mp.precision = 70
mp.mp.dps = 70


def main():
    #args = get_args()
    #if args.polfile: 
    #otherwise generate random pol with given deg, etc.
    return

if __name__ == '__main__':
    infile = sys.argv[1]

    with open(infile) as f:
        L = f.read().split('\n')

    while L[0][0] == '!': print(L.pop(0))

    L = [l for l in L if l.strip() != ''] 
    pol_type = L.pop(0)

    #case dri 
    if pol_type == 'dri':
        L.pop(0) # 0 line. 
        deg = int(L.pop(0))
    
        L.reverse()
        p_coeffs = [float(l) for l in L]
        p_degs = [i for i in range(deg+1)].reverse()
        p = lambda x: mp.polyval(p_coeffs, x)
        dp_coeffs = [p_coeffs[i]*(deg-i) for i in range(len(p_coeffs)-1)]
        dp = lambda x: mp.polyval(dp_coeffs, x)
        dp_degs = [d-1 for d in p_degs] 
        if -1 in dp_degs: 
            loc = dp_degs.index(-1) 
            dp_degs.pop(loc)
    
    #case sri
    if pol_type == 'sri':
        L.pop(0) # 0 line. 
        deg = int(L.pop(0))
        num_terms = int(L.pop(0))
        p_degs = list(map(float, L[::2]))
        p_coeffs = list(map(float, L[1::2]))
        p = lambda x: mp.fsum(list(map( lambda y,z: mp.fmul(y, mp.power(x,z)), p_coeffs, p_degs)))
        dp_degs = [d-1 for d in p_degs] 
        dp_coeffs = [p_coeffs[i]*p_degs[i] for i in range(len(p_degs))]

        #get rid of the 0-ed out constant term
        #not strictly required, since the coeff should be 0 for the term that disappears
        if -1 in dp_degs: 
            loc = dp_degs.index(-1) 
            dp_degs.pop(loc)
            dp_coeffs.pop(loc)
        dp = lambda x: mp.fsum(list(map( lambda y,z: mp.fmul(y, mp.power(x,z)), dp_coeffs, dp_degs)))

    cmd_str = 'mpsolve '+infile
    Cmd = shlex.split(cmd_str)
    output = subprocess.check_output(Cmd).strip().split(b'\n')
    num_pat = b'\((.+), (.+)\)'
    roots = []
    radii = []
    for i in range(deg):
        r = re.match(num_pat, output[i]).groups()
        rc = mp.mpc(real=float(r[0]), imag=float(r[1]))
        roots.append(rc)
        rad = mp.fabs(rc) 
        radii.append(rad)

    r1 = max(radii)
    rd = min(radii)
    print(p_coeffs)
    print(p_degs)
    print(dp_coeffs)
    print(dp_degs)
    main()
