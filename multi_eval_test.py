import cmath 
import random 
import numpy as np
from multi_eval.multi_eval import mult_eval 

#global parameters are set here
random.seed()
low    = -2**1
high   =  2**1
degree =  2**2

#Generates random complex numbers in the range low,high
def gen_c_num():
	a = random.uniform(low, high)
	b = random.uniform(low, high)
	return complex(a,b)

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
	return poly,roots

poly, roots = gen_rand_poly(degree)
der = np.polyder(poly)
print(poly)
print(der)
print(roots)
print(mult_eval(poly,[1,2,-1,0]))
poly = [1/5040, 1/180, 23/360, 7/18, 967/720, 469/180, 363/140, 1] 
der = np.polyder(poly)
print(poly)
print(der)
print(mult_eval(poly,[1,2,3,4]))
