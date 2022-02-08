from DLG_alg_mpmath import DLG_rational_form, DLG 
import numpy as np
import mpmath as mp

#complex number
j = complex(0, 1)
# define test polynomials and their derivatives
P = [lambda x:(x**2-4)*(x**2+4)]#, lambda x:x**2+2]
dP = [lambda x:4*x*(x**2)]#, lambda x: 2*x]
# define number of tests
N = 100
print("number of tests in a block = %s" % N)

for k in range(len(P)):
    p = P[k]
    dp = dP[k]

    angle = np.random.uniform(low=0, high=1)
    x = np.exp(angle*2*np.pi*j)
    for l in range(1,13):   
        print("l = %s" % l)
        # loop and record the max relative difference
        max_abs_rat = 0
        max_rel_rat = 0
        for d in range(2,16):
            max_delta = 7*10**(-d)
            print("max_delta = %s" % max_delta)
            for n in range(2,N):
                # pick a random point on the unit circle
                pert_r = np.random.uniform(low=-max_delta, high=max_delta)
                pert_i = np.random.uniform(low=-max_delta, high=max_delta)
                x_pert = x + pert_r + pert_i*j
                delta_x = np.absolute(pert_r + pert_i*j)
                if delta_x == 0: print("ERROR: zero perturbation")
    
                f = DLG(p, dp, x, l)
                f_pert = DLG(p, dp, x_pert, l)
                delta_f = np.absolute(f_pert - f)
                abs_rat = delta_f*np.reciprocal(delta_x)
                max_abs_rat = max(max_abs_rat, abs_rat)
                if f == 0:
                    #print("ERROR: zero outcome")
                    continue
                rel_rat = delta_f*np.reciprocal(np.absolute(f)*delta_x)*np.absolute(x)
                max_rel_rat = max(max_rel_rat, rel_rat)
            
            # aggregate the data
            print("estimated absolute cond num: %s" % max_abs_rat)
            print("estimated relative cond num: %s" % max_rel_rat)
        print("")
