from random_poly_gen import *
import math
import time

log_d  =  6
d      =  2**log_d

p, dp, p_rev, dp_rev, coeffs = gen_rand_poly(d)
roots = mp.polyroots(coeffs)
min = math.inf
max = -math.inf
for i in range(len(roots)):
	if float(mp.fabs(roots[i])) < min:
		min = i
	if float(mp.fabs(roots[i])) > max:
		max = i

rt     = roots[min]
rev_rt = roots[min]

#l_max is the maximum l we will try and log_epsilon_max is the -log_2 of the
#smallest epsilon to zero that we will try
l	= int(math.log2(d))


#x is the point which defines a line to 0 on which we are taking a limit
angle = mp.rand()
x = mp.expjpi(angle*2)

start = time.time()

extra_precision = int(l/3)
mp.mp.dps 		= precision + 2**extra_precision
e 				= int((precision + extra_precision)/2)


		#print(mp.mp)
print("l=%s, e=%s" % (l,e))
approx = div(d,mp.fabs(DLG(p,dp,sub(0,mul(x,mp.power(2,-e))),l)))
print("approx=",mp.fabs(approx))
real = mp.power(rt,mp.power(2,l))
print("radius=", mp.fabs(real))
print("error=", mp.fabs((mp.fabs(real))-(mp.fabs(approx)))/(mp.fabs(real)))
print("approx_root=",mp.root(mp.fabs(approx), mp.power(2,l)))
print("radius_root=", mp.fabs(rt))
print("error_root=", mp.fabs(sub(mp.fabs(mp.root(mp.fabs(approx), mp.power(2,l))), mp.fabs(rt))))


print("l=%s, e=%s" % (l,e))
approx = div(DLG(p_rev,dp_rev,sub(0,mul(x,mp.power(2,-e))),l),d)
print("approx=",mp.fabs(approx))
real = mp.power(rt,mp.power(2,l))
print("radius=", mp.fabs(real))
print("error=", mp.fabs((mp.fabs(real))-(mp.fabs(approx)))/(mp.fabs(real)))
print("approx_root=",mp.root(mp.fabs(approx), mp.power(2,l)))
print("radius_root=", mp.fabs(rt))
print("error_root=", mp.fabs(sub(mp.fabs(mp.root(mp.fabs(approx), mp.power(2,l))), mp.fabs(rt))))
