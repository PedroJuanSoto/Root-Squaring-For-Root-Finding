import os
import sys
import shlex
import subprocess
import re
import argparse
from DLG_alg_mpmath import DLG_rational_form, DLG, add, sub, mul, div, precision, mpc, mpf, pod
import mpmath as mp
import math
import time

#mp.precision = 30

# reads in given arguments

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=str, help="polynomial file in mpsolve input format")
    args = parser.parse_args()
    return args


#Should handle the following formats:
#dri drq drf dci dcq sci sri srf
def get_pols(pol_file):
    with open(pol_file) as f:
        L = f.read().split('\n')

    L = [l for l in L if l.strip() != '']
    while L[0][0] == '!': L.pop(0)
    pol_type = L.pop(0)
    L.pop(0) # 0 line. 
    deg = int(L.pop(0))

    if pol_type[2] == 'i': 
        num_lines = 1
        num_fcn = lambda x: mp.mpf(x)
    elif pol_type[2] == 'f':
        num_lines = 1
        num_fcn = lambda x: mp.mpf(x)
    elif pol_type[2] == 'q':
        num_lines = 2
        num_fcn = lambda x: mp.mpf(x)

    if pol_type[1] == 'r':
        num_comps = 1
    elif pol_type[1] == 'c':
        num_comps = 2

    #if pol_type[0] == 'd':
    #    p_degs = [i for i in range(degs + 1)]
    if pol_type[0] == 's':
        num_terms = int(L.pop(0))
        p_degs = []
        for i in range(num_terms):
            p_degs.append(int(L.pop(i*(num_lines*num_comps+1)-i)))

    #L left over is p_coeffs
    p_coeffs = list(map(num_fcn, L))
    if pol_type[2] == 'q':
        p_coeffs = [ div(x,y) for (x,y) in zip(p_coeffs[0::2], p_coeffs[1::2])]
    
    if pol_type[1] == 'c':
        p_coeffs = [ mp.mpc(real=x, imag=y) for (x,y) in zip(p_coeffs[0::2], p_coeffs[1::2])]
    
    if pol_type[0] == 'd':
        p_coeffs.reverse()
        p_degs = [deg-i for i in range(deg+1)]
        p = lambda x: mp.polyval(p_coeffs, x)
        dp_coeffs = [p_coeffs[i]*(deg-i) for i in range(len(p_coeffs)-1)]
        dp = lambda x: mp.polyval(dp_coeffs, x)

        #not strictly needed but for posterity
        dp_degs = [d-1 for d in p_degs] 
        if -1 in dp_degs: 
            loc = dp_degs.index(-1) 
            dp_degs.pop(loc)
        q = lambda x: mp.polyval(dp_coeffs.append(0), x)
    
    elif pol_type[0] == 's':
        p = lambda x: mp.fsum(list(map( lambda y,z: mp.fmul(y, mp.power(x,z)), p_coeffs, p_degs)))
        dp_degs = [d-1 for d in p_degs] 
        dp_coeffs = [p_coeffs[i]*p_degs[i] for i in range(num_terms)]

        if -1 in dp_degs:
            loc = dp_degs.index(-1)
            dp_degs.pop(loc)
            dp_coeffs.pop(loc)
        dp = lambda x: mp.fsum(list(map( lambda y,z: mp.fmul(y, mp.power(x,z)), dp_coeffs, dp_degs)))
        q = lambda x: mp.fsum(list(map( lambda y,z: mp.fmul(y, mp.power(x,z+1)), dp_coeffs, dp_degs)))
    #get rev polys
    p_rev_degs = [deg - j for j in p_degs]
    p_rev = lambda x: mp.fsum(list(map( lambda y,z: mp.fmul(y, mp.power(x,z)), p_coeffs, p_rev_degs)))
    
    #print(p_coeffs)
    #print(p_degs)
    #print(dp_coeffs)
    #print(dp_degs)
    #print(p_coeffs)
    #print(p_rev_degs)

    dp_rev_degs = [d-1 for d in p_rev_degs]
    dp_rev_coeffs = [p_coeffs[i]*p_rev_degs[i] for i in range(len(p_rev_degs))]

    if -1 in dp_rev_degs:
        loc = dp_rev_degs.index(-1)
        dp_rev_degs.pop(loc)
        dp_rev_coeffs.pop(loc)
    dp_rev = lambda x: mp.fsum(list(map( lambda y,z: mp.fmul(y, mp.power(x,z)), dp_rev_coeffs, dp_rev_degs)))
    q_rev = lambda x: mp.fsum(list(map( lambda y,z: mp.fmul(y, mp.power(x,z+1)), dp_rev_coeffs, dp_rev_degs)))

    #print(p_coeffs)
    #print(p_degs)
    #print(dp_coeffs)
    #print(dp_degs)
    #print(p_coeffs)
    #print(p_rev_degs)
    #print(dp_rev_coeffs)
    #print(dp_rev_degs)

    return deg, p, dp, p_rev, dp_rev
    #return deg, p, q, p_rev, q_rev


# runs mpsolve on the given input file and finds roots and max/min root radii
# returns the roots and the extremal root radii
def get_roots(pol_file):
    cmd_str = 'mpsolve '+ pol_file
    Cmd = shlex.split(cmd_str)
    output = subprocess.check_output(Cmd).strip().split(b'\n')
    num_pat = b'\((.+), (.+)\)'
    roots = []
    radii = []

    for i in range(len(output)):
        r = re.match(num_pat, output[i]).groups()
        rc = mp.mpc(real=float(r[0]), imag=float(r[1]))
        #roots.append(rc)
        rad = mp.fabs(rc)
        #radii.append(rad)
        roots.append((rc, rad))
    #r_min = min(radii) #min radii
    #r_max = max(radii) #max radii
    #print("min rt radius = %s" % mp.nstr(r_min, 3))
    #print("max rt radius = %s" % mp.nstr(r_max, 3))

    roots = sorted(roots, key=lambda x: x[1])
    r_min = roots[0][0]
    r_max = roots[-1][0]
    return r_min, r_max


# runs tests using the FG algorithm
def run_tests(deg, p, dp, p_rev, dp_rev, r_min, r_max):
    l = int(math.log2(deg))

    #x is the point which defines a line to 0 on which we are taking a limit
    angle = mp.rand()
    x = mp.expjpi(angle*2)
    q = lambda y: mul(y, dp(y))
    q_rev = lambda y: mul(y, dp_rev(y))

    star_min = time.time()
    print("before",mp.mp)

    extra_precision = int(l/3)
    mp.mp.dps = precision + 2**extra_precision + 300
    e = int((precision + 2**extra_precision)/2)+600

    print("after",mp.mp)

    print("l=%s, e=%s" % (l,e))
    r1 = DLG(p,q,mul(x,mp.power(2,-e)), l, e)
    print("r1 approx=", mp.nstr(r1,3))
    min_rel_err = mp.fabs(div(sub(r_min, r1),(mp.fabs(r_min))))
    print("rel_error %=", min_rel_err*100)

    print("l=%s, e=%s" % (l,e))
    rd = DLG(p_rev,q_rev,mul(x,mp.power(2,-e)), l, e)
    max_rel_err = mp.fabs(div(sub(r_max, rd),(mp.fabs(r_max))))
    print("rd approx=", mp.nstr(rd,3))
    print("rel_error %=", max_rel_err*100)

    duration = time.time()-star_min
    print("time=", duration)

    print("%s & %s & %s & %s & %s & %s & %s & %s & $[%s, %s]$\\\\" % (deg, l, mp.mp.dps, mp.nstr(rd,3), mp.nstr(min_rel_err,3), mp.nstr(r1,3), mp.nstr(max_rel_err,3), round(duration,2), mp.nstr(r_min,3), mp.nstr(r_max,3) ))

def main():
    args = get_args()
    infile = args.infile
    deg, p, dp, p_rev, dp_rev = get_pols(infile)
    #deg, p, q, p_rev, q_rev = get_pols(infile)

    l = int(math.log2(deg))+3
    extra_precision = int(l/3)
    mp.mp.dps = precision + 2**extra_precision + 100

    r_min, r_max = get_roots(infile)
    run_tests(deg, p, dp, p_rev, dp_rev, r_min, r_max)
    #run_tests(deg, p, q, p_rev, q_rev, r_min, r_max)


if __name__ == '__main__':
    main()