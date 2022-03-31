import mpmath as mp
from math import factorial
from DLG_alg_mpmath import DLG_rational_form, DLG, add, sub, mul, div, precision, mpc, mpf, pod

#each entry contains 2D list of [coeff, [(subscript, degree)]] corresponding to each term
##this might be too many levels. We'll see.
#PDP[i] = (p'/p)^(i). Need to compute up to R_{i+1} = p^(i+1)/p, for which we need to find coeffs p_0, ..., p_{i+1}
#E.g., PDP[0] = p'/p = R_1
PDP_FORMULAE = [ \
[[1, [(1,1)]]], \
[[1, [(2,1)]], [-1, [(1,2)]]], \
[[1, [(3,1)]], [-3, [(1,1), (2,1)]], [2, [(1,3)]]], \
[[1, [(4,1)]], [-4, [(1,1), (3,1)]], [-3, [(2,2)]], [12, [(1,2),(2,1)]], [-6, [(1,4)]]], \
[[1, [(5,1)]]] \
]

# Input note:
# We already put together p, dp, p_rev, dp_rev, and deg. 


#based on Alg 39
def get_cauchy_sums(p, dp, q):
    Z = [mp.expjpi(div(2*g,q)) for g in range(q)]
    R = [div(dp(z), p(z)) for z in Z]
    S = [div(mp.fsum([mul(mp.expjpi(div(2*((h+1)*k % q), q)),R[k]) for k in range(q)]), q) for h in range(q)]

    return S

# p and dp are actually scaled p_rev and dp_rev
# S[0] is incorrect but we don't use it. 
def get_cauchy_sums_rev(p, dp, q):
    Z = [mp.expjpi(div(2*g,q)) for g in range(q)]
    R = [div(dp(z), p(z)) for z in Z]
    S = [div(mp.fsum([mul(mp.expjpi(div(2*((1-h)*k % q), -q)),R[k]) for k in range(q)]), q) for h in range(q)]

    return S

#based on eqs 38
# returns coeffs p_{d-1}, ...., p_{d-q}
# in our program, we expect p to be a reverse polynomial
def get_coeffs_via_newton_iters(S):
#assumes monic poly
    #S is the list of cauchy sums s_0, ..., s_q
    coeffs = []
    q = len(S)
    for i in range(1, q):
        s = mp.fsum(list(map(lambda x,y: x*y, coeffs[::-1], S[1:i])))
        cur_coeff = mp.fdiv(mp.fadd(s, S[i]), -i)
        coeffs.append(cur_coeff)

    return coeffs

def get_Ris(tcs, l): #tcs = i+1 trailing coefficients
    R = [0]
    for i in range(1,l+2):
        R.append(div(mul(factorial(i), tcs[i]), tcs[0]))

    return R

def eval_lth_deriv(R, l):
    #if len(R) < l+1: exit(1)
    expr = PDP_FORMULAE[l]
    deriv = mp.mpf(0)
    for term in expr:
        coeff = term[0]
        monomial = mp.fprod([pod(R[t[0]], t[1]) for t in term[1]])
        deriv = add(deriv, mul(coeff, monomial))

    return deriv
    

# based on Alg 4
# computes the q trailing coeffs of p

def alg(p, dp, p_rev, dp_rev, l):
    #get the q-1 nearly-trailing coeffs via Cauchy sums + Newton Iterations

    #get the trailing term by evaluating p at 0:
    tc = mp.mpc(p(0))

    #S = get_cauchy_sums(lambda x:div(p_rev(x), tc), lambda x: div(dp_rev(x), tc), l+2)
    S = get_cauchy_sums(lambda x:div(p_rev(x), tc), lambda x: div(dp_rev(x), tc), l+2)
    coeffs = get_coeffs_via_newton_iters(S)
    #assumes monic poly

    #coeffs[i]=p_i
    coeffs = [tc]+[mul(tc, c) for c in coeffs]

    #place holder return statement
    if l > 3: return coeffs[:l], None, None

    R = get_Ris(coeffs, l)

    D = eval_lth_deriv(R, l)

    return coeffs[:l], R, D

