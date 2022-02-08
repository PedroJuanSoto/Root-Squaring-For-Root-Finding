from DLG_alg import DLG_rational_form, DLG
import numpy as np

#p is a black box polynomial and p' is its derivative.
#x = r*e^{2*pi*i*(t/u)}, and l is the number of derivatives/depth of recursion
p  = lambda x:(x**2-4)*(x**2-4)
dp = lambda x:4*(x*(x**2))

# p_rev  = lambda x:mul(sub(1,mul(4,pod(x,2))),add(1,mul(4,pod(x,2))))
# dp_rev = lambda x:mul(-1,mul(64,mul(x,pod(x,2))))

d = 4

#l_max is the maximum l we will try and log_epsilon_max is the -log_2 of the
#smallest epsilon to zero that we will try
l_max           = 10
log_epsilon_max = 10

#x is the point which defines a line to 0 on which we are taking a limit
angle = np.random.random()
x = np.exp(angle*2*np.pi*complex(0,1))


for l in range(l_max):
    print("l=",l)
    for e in range(1,20):
        print("e=", e)
        approx = d/DLG(p,dp,0-x*(2**(-e)),l)
        print("approx=",abs(approx))
        real = 2**(2**l)
        print("radius=", abs(real))
        print("error=", abs(int(abs(real))-int(abs(approx)))/int(abs(real)))

# for l in range(l_max):
#     precision *= 2
#     mp.mp.dps = precision
#     print("l=",l)
#     for e in range(1,20):
#         print("e=", e)
#         approx = div(DLG(p_rev,dp_rev,sub(0,mul(x,mp.power(2,-e))),l),div(1,d))
#         print("approx=",float(mp.fabs(approx)))
#         real = mp.power(2,mp.power(2,l))
#         print("radius=", float(mp.fabs(real)))
#         print("error=", mp.fabs(int(mp.fabs(real))-int(mp.fabs(approx)))/int(mp.fabs(real)))
