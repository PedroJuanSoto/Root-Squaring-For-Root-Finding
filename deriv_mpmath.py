import mpmath as mp

mp.mp.dps = 30
precision = 30

add = lambda x,y : mp.fadd(x, y, exact = True)
sub = lambda x,y : mp.fsub(x, y, exact = True)
mul = lambda x,y : mp.fmul(x, y, exact = True)
div = lambda x,y : mp.fdiv(x, y) #dps = precision)
pod = lambda x,y : mp.power(x, y)
mpf = lambda x   : mp.mpf(x)    #dps = precision)
mpc = lambda x,y : mp.mpc( mpf(x), mpf(y))



#the input to this function is a black box polynomial p, its derivative p',
#a representation of a complex number x = r*e^{2*pi*i*(t/u)}, and an integer l
#so that circ_DLG(p,dp,r,t,u,l) = (1/2z) p'_i(z)/p_i(z) - p'_i(z)/p_i(z) were
#z = x^(-1/2); i.e., circ_DLG(p,dp,r,t,u,l) = "the l^th difference of p'/p at x"
def DLG_rational_form(p,dp,r,t,u,l,x,d):
	root       = eps(x,l, d)
	# print("root ", list(map(lambda x: int(abs(x*1000))/1000, root)))
	base_step  = [div(dp(r),p(r)) for r in root]
	derivs     = [base_step]
	# print("deriv", list(map(lambda x: int(abs(x*1000))/1000, derivs[0])))
	for i in range(l):
		derivs.append([])
		for j in range(2**(l-i-1)):
			derivs[i+1].append(div(sub(derivs[i][2*j], derivs[i][2*j+1]),sub( root[2*j], root[2*j+1] )    ))
		root = eps(x,l-1-i, d)
		# print("root ", i, list(map(lambda x: int(abs(x*1000))/1000, root)))
		# print("deriv", i, list(map(lambda x: int(abs(x*1000))/1000, derivs[i+1])))
	return derivs[l][0]

#the input to this function is a black box polynomial p, its derivative p',
#a complex number x, and an integer l
#so that circ_DLG(p,dp,r,t,u,l) = (1/2z) p'_i(z)/p_i(z) - p'_i(z)/p_i(z) were
#z = x^(-1/2); i.e., circ_DLG(p,dp,r,t,u,l) = "the l^th difference of p'/p at x"
def DLG(p,dp,x,l,epsilon, d):
	angle = mp.arg(x)
	u     = pod(2,epsilon)
	t     = mp.fmod(angle,pod(2,-epsilon))
	r     = mpf(mp.fabs(x))
	return DLG_rational_form(p,dp,r,t,u,l,x,d)

#This gives the evaluation points
#[x-2^(l-1)*d, x-(2^(l-1) - 1)*d ,..., x-d, x+d,....x+2^(l-1)*d, x+(2^(l-1) - 1)*d]
#for computing the derivative using finite differences
def eps(x, l, d):
	if l == 1:
		return [sub(x,d),add(x,d)]
	elif l != 0:
		left  = eps(sub(x,d), l-1,d)
		right = eps(add(x,d), l-1,d)
		left.extend(right)
		return left
	else:
		return [x]
