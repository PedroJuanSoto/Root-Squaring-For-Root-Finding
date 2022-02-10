from random_root_poly_gen import *
import math
import time

log_d  =  8
d      =  2**log_d

p, dp, p_rev, dp_rev, roots = gen_rand_poly(d)
# roots = mp.polyroots(coeffs)
print(list(map(int,list(map(mp.fabs,roots)))))
min = mpf("inf")
max = roots[0]
for i in range(len(roots)):
	# print("current min=", min, "current max=", max, "current =", roots[i], "isbigger=",  mp.fabs(roots[i]) < min, "issmaller=",  mp.fabs(roots[i]) > max,)
	if mp.fabs(roots[i]) < mp.fabs(min):
		min = roots[i]
	if mp.fabs(roots[i]) > mp.fabs(max):
		max = roots[i]

rt     = min
rev_rt = max
print("min=",int(mp.fabs(rt)),"max=",int(mp.fabs(rev_rt)))
#l_max is the maximum l we will try and log_epsilon_max is the -log_2 of the
#smallest epsilon to zero that we will try
l	= int(math.log2(d))


#x is the point which defines a line to 0 on which we are taking a limit
angle = mp.rand()
x = mp.expjpi(angle*2)

start = time.time()

print("before",mp.mp)

extra_precision = int(l/3)
# mp.mp.dps 		= precision + 2**extra_precision + 100
e 				= int((precision + 2**extra_precision)/2)+200

print("after",mp.mp)

print("l=%s, e=%s" % (l,e))
approx = div(d,DLG(p,dp,sub(0,mul(x,mp.power(2,-e))),l))
print("approx=",mp.fabs(approx))
real = mp.power(rt,mp.power(2,l))
print("radius=", mp.fabs(real))
print("error=", mp.fabs(div(sub(real,approx),(mp.fabs(real)))))
print("approx_root=",mp.root(mp.fabs(approx), mp.power(2,l)))
print("radius_root=", mp.fabs(rt))
print("error_root=", mp.fabs(sub(mp.fabs(mp.root(mp.fabs(approx), mp.power(2,l))), mp.fabs(rt))))
print("rel_error_root=", int(mul(100,div(mp.fabs(sub(mp.fabs(mp.root(mp.fabs(approx), mp.power(2,l))), mp.fabs(rt))),mp.fabs(rt)))),"%")

print("l=%s, e=%s" % (l,e))
approx = div(DLG(p_rev,dp_rev,sub(0,mul(x,mp.power(2,-e))),l),d)
print("approx=",mp.fabs(approx))
real = mp.power(rev_rt,mp.power(2,l))
print("radius=", mp.fabs(real))
print("error=", mp.fabs((mp.fabs(real))-(mp.fabs(approx)))/(mp.fabs(real)))
print("approx_root=",mp.root(mp.fabs(approx), mp.power(2,l)))
print("radius_root=", mp.fabs(rev_rt))
print("error_root=", mp.fabs(sub(mp.fabs(mp.root(mp.fabs(approx), mp.power(2,l))), mp.fabs(rev_rt))))
print("rel_error_root=", int(mul(100,div(mp.fabs(sub(mp.fabs(mp.root(mp.fabs(approx), mp.power(2,l))), mp.fabs(rev_rt))),mp.fabs(rev_rt)))),"%")


print("time=",time.time()-start)
