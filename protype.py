import cmath as cm


#the kth DLG iteration of x^n - a^n as a blackbox polynomial 
def nroots_squared(a,n,k):
	def polyy(x):
		return x**n - (a**n)**(2**k)  
	return polyy

#this evaluates a polynomial p = (x-r_1)*...*(x-r_d) at x 
def p_val(p, x):
	out = 1
	for r in p:
		out *= x-r
	return out

#this performs DLG on a polynomial p = (x-r_1)*...*(x-r_d)  
def p_square(p):
	q = []
	for r in p:
		q.append(r*r)
	return q

#this returns the p_l as an array of blackbox polynomials
def roots_squared(r, k):
	p = []
	for i in r:
		p.append(i)
	for i in range(k):	
		p = p_square(p)
	def polyy(x):
		return p_val(p,x) 
	return polyy


#this is the algorithm that retrieves the root x from r = x^(2^k)
#the algorithm takes in an array p = [p_0,p_1,...,p_k] of DLG iterates
def retrieve_roots(p, r, k):
	if k == 1:
		return r
	left  =  cm.sqrt(r)
	right = -cm.sqrt(r)
	if abs(p[k](left)) >= abs(p[k](right)):
		return retrieve_roots(p, right, k-1)
	else:	
		return retrieve_roots(p, left, k-1)
		
#a test for n_roots
a = 2 	
k = 4
n = 3
r = (a+0.001)**k
p = []

for i in range(k):
	p.append(nroots_squared(a,n,k))

p.reverse

print(retrieve_roots(p,r,k-1))
 
#a test for a general polynomial
roots = [1,2,3,4]
k = 4
r = (a+0.1)**k
p = []

for i in range(k):
	p.append(roots_squared(roots,k))

p.reverse
	
print(retrieve_roots(p,r,k-1))
	 
