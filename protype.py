import cmath as cm
import math
from functools import reduce
from copy import deepcopy 

#the kth DLG iteration of x^n - a^n as a blackbox polynomial 
def nroots_squared(a,n,k):
	def polyy(x):
		return x**n - (a**n)**(2**k)  
	return polyy

#this evaluates a polynomial p = (x-r_1)*...*(x-r_d) at x 
#p is represented as p = [r_1,...,r_d]
def p_val(p, x):
	#print("i am evaluating", p, "at", x)
	return reduce(lambda x,y: x*y, list(map(lambda r: x-r, p)), 1)

#this performs DLG on a polynomial p = (x-r_1)*...*(x-r_d)  
def p_square(p):
	#print("i squared", p, "to get", list(map(lambda r: r*r, p)))
	return list(map(lambda r: r*r, p))

#this returns the p_l as an array of blackbox polynomials
#suspicious 
def roots_squared(roots, k):
	blackboxes = []
	p = deepcopy(roots) 
	ps = [p]
	blackboxes.append(lambda x: p_val(p,x))
	#print("the", 0, "th blackbox is", ps[0])
	for i in range(1,k):	
		ps.append(p_square(ps[i-1]))
		#print("the", i, "th blackbox is", ps[i])
		blackboxes.append(lambda x: p_val(ps[i],x))
	return blackboxes 


#this is the algorithm that retrieves the root x from r = x^(2^k)
#the algorithm takes in an array p = [p_0,p_1,...,p_k] of DLG iterates
def retrieve_roots(p, r, k):
	if k == 0:
		left  =  cm.sqrt(r)
		right = -cm.sqrt(r)
		#print(k, p[k], left, right,  math.log(1+abs(p[k](left))), math.log(1+abs(p[k](right))))
		if abs(p[k](left)) <= abs(p[k](right)):
			#print("at k=", k, "i went left")
			return left 
		else:	
			#print("at k=", k, "i went right")
			return right 
	left  =  cm.sqrt(r)
	right = -cm.sqrt(r)
	#print(k, p[k], left, right,  math.log(1+abs(p[k](left))), math.log(1+abs(p[k](right))))
	if abs(p[k](left)) <= abs(p[k](right)):
		#print("at k=", k, "i went left")
		return retrieve_roots(p, left, k-1)
	else:	
		#print("at k=", k, "i went right")
		return retrieve_roots(p, right, k-1)

print("\n\n\n======================================\n\n\n") 
		
#a test for n_roots
a = 2 	
k = 4
n = 3
r = (a+0.001)**(2**k)
p = []

for i in range(k):
	p.append(nroots_squared(a,n,k))

#p.reverse

print(retrieve_roots(p,r,k-1))

print("\n\n\n======================================\n\n\n") 

#a test for a general polynomial
j = complex(0,1)
r1 = 1 + 0.1*j
r2 =  -2 -j
roots = [r1, r2, 3+j, 2.5+j, 2.4-1.2*j]
a = -2.5-j
e = 0.2 + 0.1*j
k = 3
r = (a + e)**(2**k)
p = roots_squared(roots,k)
fs = retrieve_roots(p,r,k-1)

print("x**2**k = ", r)
print("roots of poly =", roots)
print("root guess = ", a+e)
print("distance of solution to roots", list(map(lambda x: abs(x - fs), roots)))
print("distance of squared guess to squared roots", list(map(lambda x: abs(x**(2**k)-r), roots)))
print("final solution", fs)
	 
print("\n\n\n======================================\n\n\n") 

#a working version 
#j = complex(0,1)
#r1 = 1 + 0.1*j
#r2 =  -2 -j
#p0 = [r1, r2]
#a = 2
#k = 2
#r = (a + 1j)**(2**k)
#p1 = [r1*r1, r2*r2]
#l0 = lambda x : (x-p0[0])*(x-p0[1]) 
#l1 = lambda x : (x-p1[0])*(x-p1[1]) 
#
#p = [l0,l1]
#
##for i in range(k):
##	p.append(roots_squared(roots,k))
#print("x**2**k = ", r)
#print("final solution", retrieve_roots(p,r,k-1))
#	 
#print("\n\n\n======================================\n\n\n") 
