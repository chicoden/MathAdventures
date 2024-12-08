from tschirnhaus import *
refactored = tschirnhaus(4, 2)

def texify_term(term, coeffs):
    frags = []
    for power, coeff in zip(reversed(term[1]), coeffs):
        if power == 1:
            frags.append(coeff)

        elif power > 1:
            frags.append(coeff + "^{" + str(power) + "}")

    if len(frags) == 0:
        return str(term[0])

    joined = " ".join(frags)
    if term[0] == 1:
        return joined

    if term[0] == -1:
        return "-" + joined

    return str(term[0]) + " " + joined

def texify_sum(summation, poly_coeffs):
    result = texify_term(summation[0], poly_coeffs)
    if len(summation) == 1:
        return result

    for term in summation[1:]:
        tex = texify_term(term, poly_coeffs)
        if tex[:1] == "-":
            result += " - " + tex[1:]

        else:
            result += " + " + tex

    return result

def texify_sum_times_term(summation, term, poly_coeffs, trans_coeffs):
    sum_string = texify_sum(summation, poly_coeffs)
    term_string = texify_term([1, term], trans_coeffs)
    if term_string == "1":
        return sum_string

    if sum_string == "1":
        return term_string

    if sum_string == "-1":
        return "-" + term_string

    try:
        int(sum_string)
        return sum_string + " " + term_string

    except ValueError:
        pass

    if len(summation) == 1:
        return sum_string + " " + term_string

    return "(" + sum_string + ") " + term_string

def texify_coeff(trans_coeff, poly_coeffs, trans_coeffs):
    result = texify_sum_times_term(*trans_coeff[0], poly_coeffs, trans_coeffs)
    if len(trans_coeff) == 0:
        return result

    for summation, term in trans_coeff[1:]:
        tex = texify_sum_times_term(summation, term, poly_coeffs, trans_coeffs)
        if tex[:1] == "-":
            result += " - " + tex[1:]

        else:
            result += " + " + tex

    return result

def tschir2latex(trans, poly_coeffs, trans_coeffs):
    result = ""
    for trans_coeff, name in zip(reversed(trans[:-1]), poly_coeffs):
        result += name + "' = " + texify_coeff(trans_coeff, poly_coeffs, trans_coeffs) + "\n\n"

    return result[:-2]
