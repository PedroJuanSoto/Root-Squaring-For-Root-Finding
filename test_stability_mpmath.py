from DLG_alg_mpmath import DLG_rational_form, DLG, add, sub, mul, div, precision, mpc, mpf
import mpmath as mp 

#complex number
j = mpc(0, 1)
# define test polynomials and their derivatives
P = [lambda x:mul(sub(mp.power(x,2),4),add(mp.power(x,2),4))]#, lambda x:x**2+2]
dP = [lambda x:mul(4,mul(x,mp.power(x,2)))]#, lambda x: 2*x]
# define number of tests
N = 100
print("number of tests in a block = %s" % N)
gen_rand = lambda low, high : mpf((high - low)*mp.rand() + low) 

for k in range(len(P)):
	p = P[k]
	dp = dP[k]

	angle = mp.rand()
	x = mp.expjpi(angle*2)
	for l in range(1,13):   
		print("l = %s" % l)
		# loop and record the max relative difference
		max_abs_rat = 0
		max_rel_rat = 0
		for d in range(2,16):
			precision += 2
			mp.dps = precision 
			print(mp.mp)
			max_delta = 7*10**(-d)
			print("max_delta = %s" % max_delta)
			for n in range(2,N):
				# pick a random point on the unit circle
				pert_r = mpf(gen_rand(low=-max_delta, high=max_delta))
				pert_i = mpf(gen_rand(low=-max_delta, high=max_delta))
				x_pert = add(x, add(pert_r, mul(pert_i,j)))
				delta_x = mpf(mp.fabs(add(pert_r,  pert_i*j)))
				if delta_x == 0: print("ERROR: zero perturbation")
	
				f = DLG(p, dp, mpc(x.real, x.imag), l)
				f_pert = DLG(p, dp, mpc(x_pert.real, x_pert.imag), l)
				delta_f = mpf(mp.fabs(sub(f_pert,  f)))
				abs_rat = div(delta_f, delta_x)
				max_abs_rat = max(max_abs_rat, abs_rat)
				if f == 0:
					#print("ERROR: zero outcome")
					continue
				rel_rat = div(mul(delta_f,mpf(mp.fabs(x))),mul(mpf(mp.fabs(f)),delta_x)),
				max_rel_rat = max(max_rel_rat, rel_rat[0])
			
			# aggregate the data
			print("estimated absolute cond num: %s" % max_abs_rat)
			print("estimated relative cond num: %s" % max_rel_rat)
		print("")
