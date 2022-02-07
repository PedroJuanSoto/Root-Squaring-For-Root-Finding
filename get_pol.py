#e.g. python3 get_pol.py mpsolve_pol/src/tests/unisolve/chrma22.pol
import os
import sys
import mpmath as mp

if __name__ == '__main__':
    infile = sys.argv[1]
    with open(infile) as f:
        L = f.read().split('\n')

    while L[0][0] == '!': L.pop(0)

    #case dri 
    L = [l for l in L if l.strip() != ''] 
    pol_type = L.pop(0)
    if pol_type == 'dri':
        L.pop(0) # 0 line. 
        deg = int(L.pop(0))
        #print("deg: ", deg)
    
        L.reverse()
        p_coeffs = [int(l) for l in L]
        #print(p_coeffs)
        p = lambda x: mp.polyval(p_coeffs, x)
        dp_coeffs = [p_coeffs[i]*(deg-i) for i in range(len(p_coeffs)-1)]
        #print(dp_coeffs)
        dp = lambda x: mp.polyval(dp_coeffs, x)
    
    #Next: case sri #sparse
