from math import ceil, log
from numpy.fft import fft, ifft
import numpy as np
from copy import deepcopy

#FFT polynomial multiplication
def poly_mult(p1, p2):
	if not p1 or not p2:
		return []
	deg1 = len(p1) - 1 
	deg2 = len(p2) - 1
	d = deg1 + deg2 + 1
	U = fft(poly_make_deg_eq(p1, d)[::-1])
	V = fft(poly_make_deg_eq(p2, d)[::-1])
	res = list(ifft(U*V).real)
	return [round(x,14) for x in res[::-1]]

#multiply a polynomial by a scalar
def poly_scalar_mult(a, p):
	return [a*pi for pi in p]

#subtract two polynomials
def poly_sub(u, v):
	d = max(len(u), len(v))
	return poly_stand_form([a-b for a,b in zip(poly_make_deg_eq(u, d), poly_make_deg_eq(v, d))])

#add two polynomials
def poly_add(u, v):
	d = max(len(u), len(v))
	return poly_stand_form([a+b for a,b in zip(poly_make_deg_eq(u, d), poly_make_deg_eq(v, d))])

#pad the polynomial p with zeros so that it has new degree = d
def poly_make_deg_eq(p, d):
	return [0] * (d-len(p)) + list(p)

#gets rid of the leading zeroes of p
def poly_stand_form(p):
	for i,a in enumerate(p):
		if a != 0:
			return p[i:]
	if 0 in p: 
		return [0]
	else:
		return []	

#fast polynomial division using newton raphson method
def poly_div(p, q):
	a = poly_deg(p)
	b = poly_deg(q)
	if a < b:
		return [0], p
	m = a-b
	y = newt_raph_inv(poly_rev(q,b),m+1)
	x = poly_rev(p,a)
	f = mod_by_x_to_the_i(poly_mult(x,y),m+1)	
	g = poly_rev(f,m)
	r = poly_sub(p,poly_mult(g,q))	
	return g,r

#newton raphson method for finding f^-1 mod x^i
def newt_raph_inv(p,i):
	g = [1]
	r = 1
	while 2**r < i:
		r += 1 
	for j in range(1,r+1):
		l = poly_scalar_mult(2,g)
		u = poly_mult(g,g)
		m = poly_mult(p,u)
		g = mod_by_x_to_the_i(poly_sub(l,m),2**j)
	return mod_by_x_to_the_i(g, i)

#returns p mod x^i for a polynomial p
def mod_by_x_to_the_i(p,i):
	return p[-i:]	

#revereses the coefficients of a polynomial
def poly_rev(p,k):
	#return rat_to_poly(rev(poly_to_rat(p),k))
	if len(p) == k:
		rev_p = deepcopy(p) 
		rev_p.reverse()
		return rev_p
	return poly_rev(poly_stand_form(p),k)

#returns the degree of a polynomial
def poly_deg(p):
	n = len(p)
	for i in range(n):
		if p[n-1-i] != 0:
			return n-1-i
	return 0

