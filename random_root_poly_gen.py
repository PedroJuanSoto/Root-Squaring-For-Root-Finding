from DLG_alg_mpmath import DLG_rational_form, DLG, add, sub, mul, div, precision, mpc, mpf, pod
import mpmath as mp
from copy import deepcopy

mp.mp.dps = 300



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

#Generates random complex numbers in the range low,high
def gen_c_num():
	a = gen_rand(low, high)
	b = gen_rand(low, high)
	return complex(a,b)

#Generates a random ploynomial of degree d whose roots are generated uniformly randomly
def gen_rand_poly(d):
	roots = []
	for i in range(d):
		roots.append(gen_c_num())
	linear_factors = []
	for r in roots:
		linear_factors.append([1,-r])
	poly = linear_factors.pop()
	for p in linear_factors:
		poly = mul_poly_lin_factor(poly,p)
	coeffs = poly
	coeffs_rev = deepcopy(coeffs)
	coeffs.reverse()
	def p_rev(x):
		return eval(coeffs,x)
	def dp_rev(x):
		return eval_deriv(coeffs,x)
	def p(x):
		return eval(coeffs_rev,x)
	def dp(x):
		return eval_deriv(coeffs_rev,x)
	return p, dp, p_rev, dp_rev, roots


def mul_poly_lin_factor(poly,lin):
	left = deepcopy(poly)
	left.append(0)
	out = [poly[0]]
	for i in range(len(poly)):
		out.append(add(left[i+1],mul(lin[1],poly[i])))
	return out
