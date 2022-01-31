from DLG_alg_mpmath import DLG_rational_form, DLG, add, sub, mul, div, precision, mpc, mpf
import mpmath as mp 
import matplotlib.pyplot as plt


#complex number
j = mpc(0, 1)
# define test polynomials and their derivatives
P = [lambda x:mul(sub(mp.power(x,2),4),add(mp.power(x,2),4))]#, lambda x:x**2+2]
dP = [lambda x:mul(4,mul(x,mp.power(x,2)))]#, lambda x: 2*x]
# define number of tests
N = 100

r = 2

min_x = 0
max_x = 0.25
num_x = 3
X = mp.linspace(min_x, max_x, num_x) #key points for evaluation, re: roots of p

min_l = 1
max_l = 12
L = [i for i in range(min_l, max_l+1, 3)]

min_d = 15
max_d = 40
D = [i for i in range(min_d, max_d+1)]

print("number of tests in a block = %s" % N)
gen_rand = lambda low, high : mpf((high - low)*mp.rand() + low) 

for k in range(len(P)):
	p = P[k]
	dp = dP[k]
	abs_res = []
	rel_res = []

	for i in range(num_x):
		abs_res.append([])
		rel_res.append([])
		angle = X[i]#mp.rand()
		x = mul(mp.expjpi(angle*2),r)
		#for l in range(1,13):   
		for l in L:
			abs_res[-1].append([])
			rel_res[-1].append([])
			print("l = %s" % l)
			# loop and record the max relative difference
			max_abs_rat = 0
			max_rel_rat = 0
			#for d in range(0,41):
			f = DLG(p, dp, mpc(x.real, x.imag), l)
			for d in D:
				precision += 2
				mp.dps = precision 
				print(mp.mp)
				max_delta = 0.7*(2**(-d))
				#print("max_delta = %s" % max_delta)
				print("x = %s, l = %s, d = %s" % (x, l, d))

				for n in range(2,N):
					# pick a random point on the unit circle
					pert_r = mpf(gen_rand(low=-max_delta, high=max_delta))
					pert_i = mpf(gen_rand(low=-max_delta, high=max_delta))
					x_pert = add(x, add(pert_r, mul(pert_i,j)))
					delta_x = mpf(mp.fabs(add(pert_r,  pert_i*j)))
					if delta_x == 0: print("ERROR: zero perturbation")
	
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
				abs_res[-1][-1].append(max_abs_rat)
				rel_res[-1][-1].append(max_rel_rat)
				print("estimated absolute cond num: %s" % max_abs_rat)
				print("estimated relative cond num: %s" % max_rel_rat)
			print("")

fig, axs = plt.subplots(nrows=num_x, ncols=1, figsize=(5,5)) #constrained_layout=True, figsize=(2,3))
for ax in axs:
	ax.remove()
gridspec = axs[0].get_subplotspec().get_gridspec()
subfigs = [fig.add_subfigure(gs) for gs in gridspec]

for i in range(num_x):
	subfig = subfigs[i]
	subfig.suptitle(f'x = {mul(mp.expjpi(angle*2),r)}') #'angle = {X[i]}')

	axs_s = subfig.subplots(nrows=1, ncols=2)
	#axs_s[0].set_title(f'absolute')
	axs_s[0].set_xlabel('-ln(error bound)')
	for j in range(len(L)):
		axs_s[0].plot(abs_res[i][j])

	#axs_s[1].set_title(f'relative')
	axs_s[1].set_xlabel('-ln(error bound)')
	for j in range(len(L)):
		axs_s[1].plot(rel_res[i][j])

plt.show()
