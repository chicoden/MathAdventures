# TODO: solve for complex roots too (maybe have seperate, only real roots solvers?)
import numpy as np

from .arithmetic import *
from .calculus import *
from .analysis import *

# x_n+1 = x_n - P(x_n) / P'(x_n)
def newtons_method(poly, start, iters, max_precis=None):
    x = start
    deriv = diff_poly(poly)
    for step in range(0, iters):
        y = eval_poly(poly, x)
        slope = eval_poly(deriv, x)

        if slope == 0: break
        x -= y / slope

        if max_precis is not None:
            if abs(y) < max_precis:
                break

    return x

# x_n+1 = x_n - 2P(x_n)P'(x_n) / (2P'(x_n)^2 - P(x_n)P''(x_n))
def halleys_method(poly, start, iters, max_precis=None):
    x = start
    deriv = diff_poly(poly)
    deriv2 = diff_poly(deriv)
    for step in range(0, iters):
        y = eval_poly(poly, x)
        slope = eval_poly(deriv, x)
        curve = eval_poly(deriv2, x)
        denom = 2 * slope * slope - y * curve

        if denom == 0: break
        x -= 2 * y * slope / denom

        if max_precis is not None:
            if abs(y) < max_precis:
                break

    return x

# Almost surefire
# G = P'(x_n) / P(x_n)
# H = G^2 - P''(x_n) / P(x_n)
# x_n+1 = x_n - deg/(G +/- sqrt((deg - 1)(deg * H - G^2)))
def laguerres_method(poly, start, iters, max_precis=None):
    x = start
    n = len(poly) - 1
    deriv = diff_poly(poly)
    deriv2 = diff_poly(deriv)
    for step in range(0, iters):
        y = eval_poly(poly, x)
        G = eval_poly(deriv, x) / y
        H = G * G - eval_poly(deriv2, x) / y
        denom = G + np.sign(G) * np.sqrt((n - 1) * (n * H - G * G))

        if denom == 0: break
        x -= n / denom

        if max_precis is not None:
            if abs(y) < max_precis:
                break

    return x

# x1_n+1 = x1_n - f(x1_n) / (x1_n - x2_n)(x1_n - x3_n)...
# x2_n+1 = x2_n - f(x2_n) / (x2_n - x1_n)(x2_n - x3_n)...
# ...
def durand_kerner_method(poly, start, iters, max_precis=None):
    x = start
    for step in range(0, iters):
        for idx, root in enumerate(x):
            denom = 1
            for other_idx, other_root in enumerate(x):
                if other_idx != idx:
                    denom *= root - other_root

            if denom == 0: break
            x[idx] -= eval_poly(poly, root) / denom

        if max_precis is not None:
            max_err = 0
            for root in x:
                max_err = max(max_err, abs(eval_poly(poly, root)))

            if max_err < max_precis:
                break

    return x

def bisection_method(p, left, right, iters, max_precis=None):
    for step in range(0, iters):
        guess = (left + right) / 2
        guess_val = eval_poly(p, guess)
        left_val = eval_poly(p, left)
        if np.sign(guess_val) == -np.sign(left_val):
            right = guess

        else:
            left = guess

        if max_precis is not None:
            if abs(guess_val) < max_precis:
                break

    return (left + right) / 2

def solve_poly(p, iters, min_bound=None, max_bound=None, max_precis=None):
    if min_bound is None or max_bound is None:
        # TODO: improve bound
        bound_radius = cauchy_root_bounding_disc(p)
        if min_bound is None: min_bound = -bound_radius
        if max_bound is None: max_bound = bound_radius

    # Map root in [x_min, x_max] interval to the [0, 1] interval
    p_mapped = [0]
    term = [1]
    var = [min_bound, max_bound - min_bound]
    for coeff in p:
        p_mapped = poly_add(p_mapped, [coeff * term_coeff for term_coeff in term])
        term = poly_mul(term, var)

    # Isolate roots
    intervals = [(0, 0, p_mapped)]
    isolated = []
    while len(intervals) > 0:
        start, level, q = intervals.pop(0)

        # Safety check before reciprocating the roots
        if q[0] == 0:
            q = q[1:]
            isolated.append((start, level, 0, None))

        # Make a change of variable to 1 / (x + 1)
        # This maps the roots in the [0, 1] interval to the interval [0, ∞]
        q_mapped = [0]
        term = [1]
        for coeff in reversed(q):
            q_mapped = poly_add(q_mapped, [coeff * term_coeff for term_coeff in term])
            term = poly_mul(term, [1, 1])

        max_roots = descartes_rule_of_signs(q_mapped)

        # A root has been isolated in this interval!
        if max_roots == 1:
            isolated.append((start, level, 1, q))

        # Not yet isolated, subdivide interval
        if max_roots > 1:
            level += 1 # Scale the interval
            start *= 2 # Unscale start

            # Unscale q
            coeff_scale = 1
            for deg in range(0, len(q)):
                q[-deg - 1] *= coeff_scale
                coeff_scale *= 2

            intervals.append((start, level, q))

            # Shift q
            q_new = [0]
            term = [1]
            for coeff in q:
                q_new = poly_add(q_new, [coeff * term_coeff for term_coeff in term])
                term = poly_mul(term, [1, 1])

            intervals.append((start + 1, level, q_new))

    # Refine intervals
    roots = []
    for start, level, width, q in isolated:
        if width > 0:
            for i in range(0, iters):
                if max_precis is not None:
                    error = width / (1 << level)
                    if error < max_precis:
                        break

                level += 1 # Scale the interval
                start *= 2 # Unscale start
        
                # Unscale q
                coeff_scale = 1
                for deg in range(0, len(q)):
                    q[-deg - 1] *= coeff_scale
                    coeff_scale *= 2

                if (q[0] < 0) == (sum(q) < 0):
                    # Sign is the same at both ends of the left half of the interval
                    # Therefore the root is in the right half and we need to shift q and the interval
                    q_new = [0]
                    term = [1]
                    for coeff in q:
                        q_new = poly_add(q_new, [coeff * term_coeff for term_coeff in term])
                        term = poly_mul(term, [1, 1])

                    q = q_new
                    start += 1

        roots.append(var[0] + var[1] * (start + width / 2) / (1 << level))

    return roots

# ax + b = 0
def solve_linear(a, b):
    return [-b / a]

# ax^2 + bx + c = 0
def solve_quadratic(a, b, c):
    dis = b * b - 4 * a * c
    if dis < 0:
        return []

    else:
        dis = np.sqrt(dis)
        denom = 2 * a
        return [(dis - b) / denom, (-dis - b) / denom]

# ax^3 + bx^2 + cx + d
def solve_cubic(a, b, c, d):
    b /= a; c /= a; d /= a; # Divide by leading coefficient to make it 1

    # Depress the cubic to x^3 + px + q by substituting x-b/3
    p = c - b * b / 3
    q = 2 * b * b * b / 27 - b * c / 3 + d
    offs = b / 3

    tp, hq = p / 3, q / 2
    dis = tp * tp * tp + hq * hq # Cubic discriminant
    if dis < 0: # Cardano's formula in complex numbers
        r = 2 * np.sqrt(-tp)
        a = np.arctan2(np.sqrt(-dis), -hq) / 3
        x, y = r * np.cos(a), r * np.sin(a)
        tx, ty = -0.5 * x, np.sqrt(0.75) * y
        return [x - offs, tx - ty - offs, tx + ty - offs]

    else: # Cardano's formula in real numbers
        dis = np.sqrt(dis)
        return [np.cbrt(dis - hq) + np.cbrt(-dis - hq) - offs]

# ax^4 + bx^3 + cx^2 + dx + e
def solve_quartic(a, b, c, d, e):
    b /= a; c /= a; d /= a; e /= a; # Divide by leading coefficient to make it 1

    # Depress the quartic to x^4 + px^2 + qx + r by substituting x-b/4
    bb = b * b
    p = (8 * c - 3 * bb) / 8
    q = (8 * d - 4 * c * b + bb * b) / 8
    r = (256 * e - 64 * d * b + 16 * c * bb - 3 * bb * bb) / 256

    # Solve for a root to (t^2)^3 + 2p(t^2)^2 + (p^2 - 4r)(t^2) - q^2 which resolves the
    # system of equations relating the product of two quadratics to the depressed quartic
    ra =  2 * p
    rb =  p * p - 4 * r
    rc = -q * q

    # Depress using the method above
    ru = ra / 3
    rp = rb - ra * ru
    rq = rc - (rb - 2 * ra * ra / 9) * ru

    rp3, rq2 = rp / 3, rq / 2
    rh = rp3 * rp3 * rp3 + rq2 * rq2
    if rh > 0: # Use Cardano's formula in the case of one real root
        rh = np.sqrt(rh)
        lmbda = np.cbrt(-rh - rq2) + np.cbrt(rh - rq2) - ru

    else: # Use complex arithmetic in the case of three real roots
        rm = np.sqrt(-rp / 3)
        lmbda = -2 * rm * np.sin(np.arcsin(3 * rq / (2 * rp * rm)) / 3) - ru

    # Solve two quadratics factored from the quartic using the cubic root
    roots = []
    if lmbda < 0: return roots
    else:
        t = np.sqrt(lmbda) # Because we solved for t^2 but want t
        alpha, beta = 2 * q / t, lmbda + ra

        u = b / 4
        t /= 2

        z = -alpha - beta
        if z > 0:
            z = np.sqrt(z) / 2
            h = +t - u
            roots += [h + z, h - z]

        w = +alpha - beta
        if w > 0:
            w = np.sqrt(w) / 2
            h = -t - u
            roots += [h + w, h - w]

        return roots
