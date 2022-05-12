from DLG_alg_mpmath import DLG_rational_form, DLG, add, sub, mul, div, precision, mpc, mpf, pod
import mpmath as mp
import math



c = 5

f  = lambda x : 2*mp.sin(x)
dp = lambda x : 2*mp.cos(x)
q  = dp
id = lambda x : 1


#l_max is the maximum l we will try and log_epsilon_max is the -log_2 of the
#smallest epsilon to zero that we will try
#e = approx parameter
l = 1
e = 100

#x is the point which defines a line to 0 on which we are taking a limit
angle = mp.rand()
x = mp.expjpi(angle*2)


print("mp.mp",mp.mp)
print("l=%s, e=%s" % (l,e))
approx = DLG(id, f, x*0.000000000001, l, e)
print("real_deriv   =", dp(0))
print("approx_derov =",mp.fabs(approx))
