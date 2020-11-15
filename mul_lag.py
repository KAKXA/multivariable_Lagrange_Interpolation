from sage.all import *

# necessary information:
# the size of field
# F: GF(q)
# R: polynomial ring over F
# H: the value that x can be
# h: size of H
q = 7
M = 3
F = GF(q)
R = F['x']
H = [0,1,2]
h = len(H)

def get_vector(x): 
    # develop a m-dim vector from a number between [0,h ^ m - 1]
    l = []
    for i in range(M):
        l.append(x % h)
        x //= h
    return l

def single_lag(xy_pairs):
    # single variate lagrange interpolation
    poly = R.lagrange_polynomial(xy_pairs)
    coe = poly.coefficients()
    exp = poly.exponents()
    res = [0] * h
    # turn it into the from of M-dim vector
    for i in range(len(exp)):
        res[exp[i]] = coe[i]
    return res

def mul_lag(m, points, res_coes, order):
    # m: current number of variables
    # points: (0 ~ h ** m - 1) points' y coordinate 
    # res_coes: result coeffients, which is a h^M-dim vector
    # order: tell the last recursion where to put its result coeffients into res_coes
    coes = [([None] * h) for i in range(h ** (m - 1))]
    for i in range(h ** (m - 1)):
        xy_pairs = [None] * h
        for j in range(h):
            xy_pairs[j] = (j, points[i * h + j])
        coes[i] = single_lag(xy_pairs)
    # the end of recursion
    if m == 1:
        for i in range(h):
            res_coes[order + i * (h ** (M - 1))] = coes[0][i]
    else:
        for i in range(h):
            new_points = [None] * (h ** (m - 1))
            for j in range(h ** (m - 1)):
                new_points[j] = coes[j][i]
            mul_lag(m - 1, new_points, res_coes, order + i * h ** (M - m))

def print_poly(coes):
    # turn h^M-dim coeffient vector into polynomial form
    l = len(coes)
    flag = False
    for i in range(l):
        if coes[i] != 0:
            if flag:
                print(' + ',end='')
            print(coes[i],end='')
            flag= True

            v = get_vector(i)
            for i in range(M):
                if v[i] != 0:
                    if v[i] == 1:
                        print('*x' + str(i), end='')
                    else:
                        print('*x' + str(i) + '^' + str(v[i]),end='')


if __name__ == "__main__":
    # construct a original polynomial that is the goal of lagrange interpolation
    X = PolynomialRing(F, M, 'X').gens()
    # this is what we are expected to get
    f = 1 + 2 * X[0]  ** 2 + 3 * X[1] ** 2 + 6 * X[0] * X[1] ** 2 * X[2]

    res_coes = [0] * (h ** M)
    points = [None] * (h ** M)
    # get points
    for i in range(h ** M):
        points[i] = f(*(get_vector(i)))

    mul_lag(M, points, res_coes, 0)
    print('The original polynomial:')
    print(f)
    print('The polynomial generate from lagrange interpolation:')
    print_poly(res_coes)