from DLG_alg_mpmath import DLG_rational_form, DLG, add, sub, mul, div, precision, mpc, mpf, pod
import mpmath as mp
import math


p_coeffs = [1,-10,35,-50,24]
q_coeffs =  [4,-30, 70, -50 , 0 ]

p = lambda x : mp.polyval(p_coeffs,x)
q = lambda x : mp.polyval(q_coeffs,x)


log_d  =  8
d      =  2**log_d

#l_max is the maximum l we will try and log_epsilon_max is the -log_2 of the
#smallest epsilon to zero that we will try
l	= int(math.log2(d))


#x is the point which defines a line to 0 on which we are taking a limit
angle = mp.rand()
x = mp.expjpi(angle*2)


print("mp.mp",mp.mp)

extra_precision = int(l/3)
# mp.mp.dps 		= precision + 2**extra_precision + 100
e 				= 64

print("l=%s, e=%s" % (l,e))
approx = DLG(p,q,mul(x,mp.power(2,-e)),l, e)
print("approx=",mp.fabs(approx))
