from DLG_alg_mpmath import DLG_rational_form, DLG, add, sub, mul, div, precision, mpc, mpf
import mpmath as mp

#p is a black box polynomial and p' is its derivative.
#x = r*e^{2*pi*i*(t/u)}, and l is the number of derivatives/depth of recursion
p  = lambda x:(x**2-4)*(x**2+4)
dp  = lambda x:4*x*(x**2)

#l_max is the maximum l we will try and log_epsilon_max is the -log_2 of the
#smallest epsilon to zero that we will try
l_max           = 10
log_epsilon_max = 10

#x is the point which defines a line to 0 on which we are taking a limit
angle = mp.rand()
x = mp.expjpi(angle*2)

#We want to "test" if the limit to 0 makes sense, to do this we check whether
#the sequence to zero to "cauchy-like" by checking if differences between points
# #in the sequnce get smaller
# new  = mpc(0,0)
# prev = mpc(0,0)
# for l in range(l_max):
#     precision += 2
#     mp.mp.dps = precision
#     for e in range(1,20):
#         new = DLG(p,dp,sub(0,mul(x,mp.power(2,-e))),l)
#         print(l,e,mp.fabs(sub(new,prev)))
#         prev = new

#similar to the previous test bu now we try to vary epsilon and l at the same
#time to see if it helps the limit converge
# new  = mpc(0,0)
# prev = mpc(0,0)
# for e in range(1,20):
#     precision *= 2
#     mp.mp.dps = precision
#     new = DLG(p,dp,sub(0,mul(x,mp.power(2,-e))),e)
#     print(e,int(mp.fabs(mp.log(sub(new,prev)))))
#     prev = new

# #print the absolute values as we approch 0
# for e in range(1,20):
#     precision *= 2
#     mp.mp.dps = precision
#     print(float(mp.fabs(DLG(p,dp,sub(0,mul(x,mp.power(2,-e))),e))))

# #print the absolute values as we approch 0
for l in range(l_max):
    precision += 2
    mp.mp.dps = precision
    print("l=",l)
    for e in range(1,20):
        print("e=,", e, "|DLG(2+2^(-e))|=",float(mp.fabs(DLG(p,dp,sub(2,mul(x,mp.power(2,-e))),l))))
