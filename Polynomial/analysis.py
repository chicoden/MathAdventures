import numpy as np

# Number of positive real roots is at most the number
# Of sign changes in the sequence of coefficients
# (excluding coefficients of 0)
def descartes_rule_of_signs(p):
    sign_delta = 0
    prev_sign = np.sign(p[0])
    for coeff in p[1:]:
        cur_sign = np.sign(coeff)
        if cur_sign != 0:
            if cur_sign == -prev_sign:
                sign_delta += 1

            prev_sign = cur_sign

    return sign_delta

def laguerre_real_root_bounds(p):
    n = len(p) - 1
    c, b, a = p[-3:]
    denom = n * a
    pseudo_dis = b * b - 2 * n * a * c / (n - 1)
    if pseudo_dis < 0:
        raise ValueError("this bound is undefined for polynomials with complex roots")

    offs = (n - 1) * np.sqrt(pseudo_dis)
    return ((-b - offs) / denom, (-b + offs) / denom)

def lagrange_root_bounding_disc(p):
    radius = 0
    for coeff in p[:-1]:
        radius += abs(coeff / p[-1])

    return max(1, radius)


def cauchy_root_bounding_disc(p):
    radius = 0
    for coeff in p[:-1]:
        radius = max(radius, abs(coeff / p[-1]))

    return 1 + radius
