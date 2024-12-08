# dP(x)/dx
def diff_poly(p):
    return [coeff * (deg + 1) for deg, coeff in enumerate(p[1:])]

# ∫P(x)dx
def int_poly(p):
    return [0] + [coeff / (deg + 1) for deg, coeff in enumerate(p)]
