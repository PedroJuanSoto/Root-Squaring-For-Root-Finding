import numpy as np

#the input to this function is a representation
#of a complex number x = e^{2*pi*i*(p/q)} so that
#circ_sq_root(p,q) = r,s and y = e^{2*pi*i*(r/s)}
#implies y^2 = x
def circ_sq_root(p,q):
	if p%q == 0:
		return 1,1
	else:
		return p, 2*q

#the input to this function is a representation
#of a complex number x = e^{2*pi*i*(p/q)} so that
#circ_neg(p,q) = t,u and y = e^{2*pi*i*(t/u)}
#implies -y = x
def circ_neg(p,q):
	if p%q == 0:
		return 1, 2
	else:
		return 2*p+q, 2*q

#the input to this function is a representation
#of a complex number x = e^{2*pi*i*(p/q)} and an integer l
#so that circ_roots_rational_form(p,q,l) = r_1,s_1;...;r_{2^l},s_{2^l}
#and y_i = e^{2*pi*i*(r_i/s_i)} implies y_i^{2^l} = x
#with the r_i,s_i ordered according to the DLG recursion
def circ_roots_rational_form(p,q,l):
	r, s  = circ_sq_root(p,q)
	t, u  = circ_neg(r,s)
	if l == 1:
		return [(r,s),(t,u)]
	elif l != 0:
		left  = circ_roots_rational_form(r,s,l-1)
		right = circ_roots_rational_form(t,u,l-1)
		left.extend(right)
		return left
	else:
		return [(p,q)]

#the input to this function is a representation
#of a complex number x = e^{2*pi*i*(p/q)} and an integer l
#so that circ_roots(p,q,l) = y_1,...,y_{2^l} where y_i^{2^l} = x
#with the y_i ordered according to the DLG recursion
def circ_roots(p,q,l):
	i = complex(0,1)
	roots = circ_roots_rational_form(p,q,l)
	return [np.exp(2*np.pi*i*(r/s)) for r,s in roots]

#the input to this function is a representation of a complex number
#x = r*e^{2*pi*i*(t/u)}, and an integer l
#so that roots(r,t,u,l) = y_1,...,y_{2^l} where y_i^{2^l} = x
#with the y_i ordered according to the DLG recursion
def roots(r,t,u,l):
	circ_root       = circ_roots(t,u,l)
	r_root          = np.power(r,1/(2**l))
	return [r_root*root for root in circ_root]

#the input to this function is a black box polynomial p, its derivative p',
#a representation of a complex number x = r*e^{2*pi*i*(t/u)}, and an integer l
#so that circ_DLG(p,dp,r,t,u,l) = (1/2z) p'_i(z)/p_i(z) - p'_i(z)/p_i(z) were
#z = x^(-1/2); i.e., circ_DLG(p,dp,r,t,u,l) = "the l^th difference of p'/p at x"
def DLG_rational_form(p,dp,r,t,u,l):
	root       = roots(r,t,u,l)
	base_step  = [dp(r)*np.reciprocal(p(r)) for r in root]
	derivs     = [base_step]
	for i in range(l):
		derivs.append([])
		for j in range(2**(l-i-1)):
			derivs[i+1].append((np.reciprocal(root[2*j])/2)*(derivs[i][2*j] - derivs[i][2*j+1]))
		root = roots(r,t,u,l-1-i)
	return derivs[l][0]

#the input to this function is a black box polynomial p, its derivative p',
#a complex number x, and an integer l
#so that circ_DLG(p,dp,r,t,u,l) = (1/2z) p'_i(z)/p_i(z) - p'_i(z)/p_i(z) were
#z = x^(-1/2); i.e., circ_DLG(p,dp,r,t,u,l) = "the l^th difference of p'/p at x"
def DLG(p,dp,x,l):
	angle     = np.angle(x)
	i         = complex(0,1)
	angle_deg = np.absolute(angle/(2*np.pi*i))
	t,u       = angle_deg.as_integer_ratio()
	r         = np.absolute(x)
	return DLG_rational_form(p,dp,r,t,u,l)


p  = lambda x:(x**2-4)*(x**2+4)
dp = lambda x:4*x*(x**2)
r  = 1
t  = 1
u  = 1
l  = 2
print(DLG_rational_form(p,dp,r,t,u,l))
x=1
print(DLG(p,dp,x,l))
