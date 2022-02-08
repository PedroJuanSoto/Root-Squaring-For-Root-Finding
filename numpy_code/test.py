from DLG_alg import DLG_rational_form, DLG

#p is a black box polynomial and p' is its derivative.
#x = r*e^{2*pi*i*(t/u)}, and l is the number of derivatives/depth of recursion
p  = lambda x:(x**2-4)*(x**2+4)
dp = lambda x:4*x*(x**2)
r  = 1
t  = 1
u  = 1
l  = 2
print(DLG_rational_form(p,dp,r,t,u,l))
x=1
print(DLG(p,dp,x,l))
