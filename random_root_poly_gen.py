from DLG_alg_mpmath import DLG_rational_form, DLG, add, sub, mul, div, precision, mpc, mpf, pod
import mpmath as mp
from copy import deepcopy

eval 	   = lambda coeffs,x : mp.polyval(coeffs,x)
eval_deriv = lambda coeffs,x : mp.polyval(coeffs,x,True)[1]

#global parameters are set here
log_random_root_radius =  6
low                    = -2**log_random_root_radius
high   	               =  2**log_random_root_radius
log_d                  =  2
degree                 =  2**log_d

#uniform random
gen_rand = lambda low, high : mpf((high - low)*mp.rand() + low)

import cmath
import random
import numpy as np


#Generates random complex numbers in the range low,high
def gen_c_num():
	a = random.uniform(low, high)
	b = random.uniform(low, high)
	return complex(a,b)

#converts numpy to mpmath
np_to_mpc = lambda z : mpc(z.real,z.imag)

#Generates a random ploynomial of degree d whose roots are generated uniformly randomly
def gen_rand_poly(d):
	roots = []
	for i in range(d):
		roots.append(gen_c_num())
	linear_factors = []
	for r in roots:
		linear_factors.append(np.poly1d([1,-r]))
	poly = np.poly1d([1])
	for p in linear_factors:
		poly *= p
	coeffs = list(map(np_to_mpc,list(poly)))
	coeffs_rev = deepcopy(coeffs)
	coeffs.reverse()
	def p_rev(x):
		return eval(coeffs,x)
		# result = mpc(0,0)
		# for i,c in enumerate(coeffs):
		# 	result = add(result,mul(c,pod(x,i)))
		# return result
	def dp_rev(x):
		return eval_deriv(coeffs,x)
		# result = mpc(0,0)
		# for i in range(1,len(coeffs)):
		# 	result = add(result,mul(mul(i,coeffs[i]),pod(x,i-1)))
		# return result
	def p(x):
		return eval(coeffs_rev,x)
		# result = mpc(0,0)
		# for i,c in enumerate(coeffs_rev):
		# 	result = add(result,mul(c,pod(x,i)))
		# return result
	def dp(x):
		return eval_deriv(coeffs_rev,x)
		# result = mpc(0,0)
		# for i in range(1,len(coeffs_rev)):
		# 	result = add(result,mul(mul(i,coeffs_rev[i]),pod(x,i-1)))
		# return result
	return p, dp, p_rev, dp_rev, roots
