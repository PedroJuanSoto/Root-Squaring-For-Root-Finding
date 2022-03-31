import newton_ids_alg as ni
import mpmath as mp

def get_test_pols(coeffs):
    d = len(coeffs)-1
    p = lambda x: mp.fsum([mp.fmul(coeffs[i], mp.power(x, d-i)) for i in range(len(coeffs))])
    dp = lambda x: mp.fsum([mp.fmul((d-i)*coeffs[i], mp.power(x, d-i-1)) for i in range(len(coeffs)-1)])
    coeffs_rev = coeffs[::-1]
    p_rev = lambda x: mp.fsum([mp.fmul(coeffs_rev[i], mp.power(x, d-i)) for i in range(len(coeffs))])
    dp_rev = lambda x: mp.fsum([mp.fmul((d-i)*coeffs_rev[i], mp.power(x, d-i-1)) for i in range(len(coeffs)-1)])

    return p, dp, p_rev, dp_rev

def pretty_print(mpc_list):
    for m in mpc_list:
        print(mp.nstr(m, 5))

mp.mp.dps = 70
print(mp.mp)

print("\n========================================================")
print("=======================UNIT TESTS=======================")
print("========================================================\n")
print("Test 1: newton's identities")
j = complex(0,1)
#powers of roots
#x = [[-1, j, -j], [1, -1, -1], [1, -j, j]]
#S = list(map(lambda x: sum(x), x))

#1st-3rd power sums
#print("Cauchy sums:")
#S = [-1, -1, -1]
S = [3, 0.25, 0.5625, 0.015625] #1-3 Cauchy sums
#S = [0, 0, 0, 4]
#q = len(S)
#coeffs = []
#for i in range(1,q):
#    s = mp.fsum(list(map(lambda x,y: x*y, coeffs, S[1:i])))
#    cur_coeff = mp.fdiv(mp.fadd(s, S[i]), -i)
#    coeffs.insert(0, cur_coeff)
C = [1, -0.25, -0.25, 0.0625]
print("correct output:\n%s" % C[1:])
coeffs = ni.get_coeffs_via_newton_ids(S)
print("test output:")
print(coeffs)

print("\n========================================================\n")

print("Test 2: roots inside the unit disc")
p = lambda x: x**3 - 0.25*x**2 - 0.25*x + 0.0625 #roots = 0.5, -0.5, 0.25
d = 3
C = [1, -0.25, -0.25, 0.0625]
S = [3, 0.25, 0.5625, 0.015625] #0-3 Cauchy sums
dp = lambda x: 3*x**2 - 0.5*x - 0.25
aC = ni.get_coeffs_via_newton_ids(S)
print("\ncoeffs:\n%s" % C)
print("test output: %s" % aC)


print("\n========================================================\n")

#print("Test 2: roots on the unit circle")
#p = lambda x: x**3 + x**2 + x + 1
#dp = lambda x: 3*x**2 + 2*x + 1
#print(ni.get_coeffs_via_newton_ids(p, dp, 4)) #use q=l+1?

print("Test 3: cauchy sums")
p = lambda x: x**3 - 0.25*x**2 - 0.25*x + 0.0625 #roots = 0.5, -0.5, 0.25
d = 3
dp = lambda x: 3*x**2 - 0.5*x - 0.25

print("expected output: %s" % [3, 0.25, 0.5625, 0.015625] )
S = ni.get_cauchy_sums(p, dp, d+1)
print("test output:")
pretty_print(S)

print("\n========================================================\n")

print("Test 4: integration test (near-leading coeffs of p)")
p = lambda x: x**3 - 0.25*x**2 - 0.25*x + 0.0625 #roots = 0.5, -0.5, 0.25
d = 3
C = [1, -0.25, -0.25, 0.0625]
dp = lambda x: 3*x**2 - 0.5*x - 0.25
S = ni.get_cauchy_sums(p, dp, d+1)
aC = ni.get_coeffs_via_newton_ids(S)
print("expected output:\n%s" % C)
print("test output:")
pretty_print(aC)

print("\n========================================================\n")

print("Test 5: roots of p inside unit disc: 0.5, -0.5, 0.25")
#p = lambda x: x**3 - 0.25*x**2 - 0.25*x + 0.0625 #roots = 0.5, -0.5, 0.25
#dp = lambda x: 3*x**2 - 0.5*x - 0.25
#d = 3
#p_rev = lambda x: 0.0625*(x**3) - 0.25*(x**2) - 0.25*x + 1 #roots = 2, -2, 4
#dp_rev = lambda x: 0.1875*x**2 - 0.5*x - 0.25

C = [1, -0.25, -0.25, 0.0625]
d = len(C)-1
l = d+1
p, dp, p_rev, dp_rev = get_test_pols(C)
aC, R, D = ni.alg(p, dp, p_rev, dp_rev, d, l)

#print("Should see %d trailing coeffs:" % l)
print("expected output:\n%s" % C[::-1][:l+1])
print("test output:")
pretty_print(aC)
#print("------------------------------------------")
#aC, R, D = ni.alg(p_rev, dp_rev, p, dp, l)
#print("Should see %d leading coeffs in reverse order:" % l)
#print("exact rev coeffs:\n%s" % C[::-1][-l:])
#print("approx rev coeffs:")
#pretty_print(aC)

print("\n========================================================\n")

print("Test 6: roots of p outside the unit disc: 2, -2, 4")
C = [1, -4, -4, 16] 
d = len(C)-1
l = d+1
p, dp, p_rev, dp_rev = get_test_pols(C)

aC, R, D = ni.alg(p, dp, p_rev, dp_rev, d, l)
#print("expected output:\n%s" % C)
#print("Should see %d trailing coeffs:" % l)
print("expected output:\n%s" % C[::-1][:l+1])
print("test output:")
pretty_print(aC)

#print("------------------------------------------")
#aC, R, D = ni.alg(p_rev, dp_rev, p, dp, l)
#print("Should see %d leading coeffs in reverse order:" % l)
#print("exact rev coeffs:\n%s" % C[::-1][-l:])
#print("approx rev coeffs:")
#pretty_print(aC)

print("\n========================================================\n")

print("Test 7: roots of p inside unit disc: 0.2, 0.3, 0.4, 0.5")

C = [1, -1.4, 0.71, 0.154, 0.012]
d = len(C)-1
l = d+1
p, dp, p_rev, dp_rev = get_test_pols(C)

aC, R, D = ni.alg(p, dp, p_rev, dp_rev, d, l)
#print("Should see %d trailing coeffs:" % l)
print("expected output:\n%s" % C[::-1][:l+1])
print("test output:")
pretty_print(aC)

print("\n========================================================\n")

print("Test 8: roots of p: 2, -3, 4, -5")
C = [1, 2, -25, -26, 120]
d = len(C)-1
l = d+1
p, dp, p_rev, dp_rev = get_test_pols(C)

aC, R, D = ni.alg(p, dp, p_rev, dp_rev, d, l)
#print("Should see %d trailing coeffs:" % l)
print("expected output:\n%s" % C[::-1][:l+1])
print("test output:")
pretty_print(aC)

print("\n========================================================\n")

print("Test 8a: roots of p: -2, 3, -4, 5 (reverse of 8)")
C = [1, - 2, -25, 26, 120]
d = len(C)-1
l = d+1
p, dp, p_rev, dp_rev = get_test_pols(C)

aC, R, D = ni.alg(p, dp, p_rev, dp_rev, d, l)
#print("Should see %d trailing coeffs:" % l)
print("expected output:\n%s" % C[::-1][:l+1])
print("test output:")
pretty_print(aC)

print("\n========================================================\n")

print("Test 9: roots of p: -5, -10, 15, 20")
C = [1, -20, -175, 2750, 15000]
d = len(C)-1
l = d+1
p, dp, p_rev, dp_rev = get_test_pols(C)

aC, R, D = ni.alg(p, dp, p_rev, dp_rev, d, l)
#print("Should see %d trailing coeffs:" % l)
print("expected output:\n%s" % C[::-1][:l])
print("test output:")
pretty_print(aC)

print("\n========================================================\n")

print("Test 10: roots of p: -50, -100, 150, 200")
C = [1, -200, -17500, 2750000, 150000000]
d = len(C)-1
l = d-1
p, dp, p_rev, dp_rev = get_test_pols(C)

aC, R, D = ni.alg(p, dp, p_rev, dp_rev, d, l)
#print("Should see %d trailing coeffs:" % l)
print("expected output:\n%s" % C[::-1][:l])
print("test output:")
pretty_print(aC)

print("\n========================================================\n")

print("Test 11: roots of p: -4.3, -18, 7.6, 2.3, 5.5")
C = [1, 6.9, -194.09, 315.939, 3423.46, -7441.24]
d = len(C)-1
l = 5
p, dp, p_rev, dp_rev = get_test_pols(C)

aC, R, D = ni.alg(p, dp, p_rev, dp_rev, d, l)
#print("Should see %d trailing coeffs:" % l)
print("expected output:\n%s" % C[::-1][:l])
print("test output:")
pretty_print(aC)

print("\n========================================================\n")

print("Test 12: roots of p: 2, 2.1, 2.2, 2.3, 2.4, 2.5, -5, -5.1, -5.2, -5.3, -5.4, -5.5")
C = [1,18,63.95, -525.75, -3081.05,7455.93, 52293.2, -87424.9, -429587, 854518, 1.25278e6, -4.0937e6, 2.6615e6]
d = len(C)-1
l = d+1 #3
p, dp, p_rev, dp_rev = get_test_pols(C)

aC, R, D = ni.alg(p, dp, p_rev, dp_rev, d, l)
#print("Should see %d trailing coeffs:" % l)
print("expected output:\n%s" % C[::-1][:l])
print("test output:")
pretty_print(aC)
#print("R values:")
#pretty_print(R)
#print("(p'/p)^(%d) = %s" % (l, mp.nstr(D,3)))

print("\n========================================================\n")

print("Test 13:  roots of p both inside and outside unit circle: 0.5, -0.5, 1.5, -1.5")
C = [1, 0, -2.5, 0, 0.5625]
d = len(C)-1
l = 5
p, dp, p_rev, dp_rev = get_test_pols(C)

aC, R, D = ni.alg(p, dp, p_rev, dp_rev, d, l)
#print("Should see %d trailing coeffs:" % l)
print("expected output:\n%s" % C[::-1][:l+1])
print("test output:")
pretty_print(aC)
#print("R values:")
#pretty_print(R)
#print("(p'/p)^(%d) = %s" % (l, mp.nstr(D,3)))

print("\n========================================================\n")

print("Test 14:  roots of p both inside and outside unit circle: 0.5, 1.5, -1.5")
C = [1, -0.5, -2.25, 1.125]
d = len(C)-1
l = 10
p, dp, p_rev, dp_rev = get_test_pols(C)

aC, R, D = ni.alg(p, dp, p_rev, dp_rev, d, l)
#print("Should see %d trailing coeffs:" % l)
print("expected output:\n%s" % C[::-1][:l+1])
print("test output:")
pretty_print(aC)

print("\n========================================================\n")

print("Test 15:  roots of p outside unit circle: 1.5, -1.5")
C = [1, 0, -2.25]
d = len(C)-1
l = 10
p, dp, p_rev, dp_rev = get_test_pols(C)

aC, R, D = ni.alg(p, dp, p_rev, dp_rev, d, l)
#print("Should see %d trailing coeffs:" % l)
print("expected output:\n%s" % C[::-1][:l+1])
print("test output:")
pretty_print(aC)

print("\n========================================================\n")
print("Test N: roots of p: ")
C = [1]
d = len(C)-1
l = 3
p, dp, p_rev, dp_rev = get_test_pols(C)

aC, R, D = ni.alg(p, dp, p_rev, dp_rev, d, l)
print("Should see %d trailing coeffs:" % l)
print("expected output:\n%s" % C[::-1][:l+1])
print("test output:")
pretty_print(aC)
print("R values:")
pretty_print(R)
print("(p'/p)^(%d) = %s" % (l, mp.nstr(D,3)))
