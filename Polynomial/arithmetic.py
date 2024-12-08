# [a0, a1, a2...] = a0 + a1*x + a2*x^2...
def eval_poly(p, x):
    # Horner's method
    # a0 + x*(a1 + x*(a2 + ...))
    y = p[-1]
    for coeff in reversed(p[:-1]):
        y = y * x + coeff

    return y

# A(x) + B(x)
def poly_add(a, b):
    res = a.copy()
    for deg, coeff in enumerate(b):
        if deg == len(res):
            res.append(coeff)

        else:
            res[deg] += coeff

    return res

# A(x) - B(x)
def poly_sub(a, b):
    res = a.copy()
    for deg, coeff in enumerate(b):
        if deg == len(res):
            res.append(-coeff)

        else:
            res[deg] -= coeff

    return res

# A(x) * B(x)
def poly_mul(a, b):
    res = []
    for deg_a, coeff_a in enumerate(a):
        for deg_b, coeff_b in enumerate(b):
            term_deg = deg_a + deg_b
            res_coeff = coeff_a * coeff_b
            if term_deg == len(res):
                res.append(res_coeff)

            else:
                res[term_deg] += res_coeff

    return res

# A(x) / B(x)
def poly_div(a, b):
    res = []
    rem = a.copy()
    while len(rem) > len(b) - 1:
        quot = rem[-1] / b[-1]
        quot_deg = len(rem) - len(b)
        res.insert(0, quot)
        del rem[-1]
        for deg, coeff in enumerate(b[:-1]):
            rem[deg + quot_deg] -= quot * coeff

    return res, rem

# Specialized division for eliminating a single root
# P(x) / (x - root)
def eliminate_root(p, root):
    new_poly = [p[-1]]
    for coeff in reversed(p[1:-1]):
        new_poly.insert(0, new_poly[0] * root + coeff)

    return new_poly

# P(x)^n
def poly_pow(p, n):
    res = [1]
    for mul in range(0, n):
        res = poly_mul(res, p)

    return res
