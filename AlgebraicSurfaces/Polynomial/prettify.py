def stringify_poly(p):
    string = ""
    for deg, coeff in enumerate(p):
        if coeff != 0:
            term = " - " if coeff < 0 else " + "
            coeff = abs(coeff)
            if coeff != 1 or deg == 0: term += str(coeff)

            if deg > 0:
                term += "x"
                if deg > 1: term += "^" + str(deg)

            string = term + string

    string = string[3:]
    if p[-1] < 0:
        string = "-" + string

    return string
