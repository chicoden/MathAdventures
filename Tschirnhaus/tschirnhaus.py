from itertools import combinations as combos
from more_itertools import distinct_permutations as perms

def zeros(num):
    return tuple(0 for _ in range(num))

def tag_tuple(num, pos, val=1):
    return tuple(val if idx == pos else 0 for idx in range(num))

def distrib_const(const, series):
    return [[term[0] * const, term[1]] for term in series]

def mul_power_terms(left, right):
    return [left[0] * right[0], tuple(powa + powb for powa, powb in zip(left[1], right[1]))]

def simplify_power_sum(series):
    return [term for term in series if term[0] != 0]

def accumulate_series(accum, addition):
    for term in addition:
        # Search for like term
        found_like_term = False
        for index in range(len(accum)):
            if accum[index][1] == term[1]: # Is this a like term?
                accum[index][0] += term[0] # Add onto like term
                found_like_term = True
                break

        if not found_like_term:
            accum.append(term) # Add term to accumulator if there is no like term

def reduce_esym(index, deg):
    return [-1 if index % 2 == 1 else 1, tag_tuple(deg, deg - index)]

def reduce_prime_root_sum(series, esym, poly_deg):
    assert set(series) == set(perms(series[0]))
    power_set = set(series[0])
    if power_set == {0}:
        return [[1, zeros(poly_deg)]]

    if power_set == {0, 1}:
        # Reduce to an elementary symmetric polynomial, then to a coefficient
        esym_index = series[0].count(1)
        return [reduce_esym(esym_index, poly_deg)]

    semifactors = []
    term = [*series[0]]
    while max(term) > 0:
        # Determine largest possible factor
        semifactor = poly_deg - term.count(0) # esym index
        semifactors.append(semifactor)
        for reduction in esym[semifactor - 1]:
            # Check if this term can be divided out of the leading term in the series
            if all(powa >= powb for powa, powb in zip(term, reduction)):
                term = [powa - powb for powa, powb in zip(term, reduction)]
                break

    # Multiply out the approximate factorization
    stack = [*semifactors]
    product = [*esym[stack.pop() - 1]]
    while len(stack) > 0:
        product = [tuple(powa + powb for powa, powb in zip(left, right)) for left in esym[stack.pop() - 1] for right in product]

    # Subtract the original series to get a residue
    # This SHOULD be composed of simpler root power-product series which can be reduced recursively
    for term in series:
        product.remove(term)

    # Reduce product of semifactors
    semifactz = [1, zeros(poly_deg)]
    for semifactor in semifactors:
        semifactz = mul_power_terms(semifactz, reduce_esym(semifactor, poly_deg))

    # Group residues into subseries and reduce them recursively
    # Then subtract them from the semifactorization
    result = [semifactz]
    while len(product) > 0:
        term = product.pop()
        perm_id = sorted(term)

        subseries = [term]
        index = 0
        while index < len(product):
            if sorted(product[index]) == perm_id:
                subseries.append(product.pop(index))

            else:
                index += 1

        multiplicity = subseries.count(term)
        minimal_subseries = [*{*subseries}] # Remove duplicates
        assert all(subseries.count(other) == multiplicity for other in minimal_subseries)
        accumulate_series(result, distrib_const(-multiplicity, reduce_prime_root_sum(minimal_subseries, esym, poly_deg)))

    return result

def reduce_root_sum(series, esym, poly_deg):
    # Decompose the series into simpler subseries and reduce them
    reduced = []
    while len(series) > 0:
        term = series.pop()
        perm_id = sorted(term)

        # Collect terms with the same multiset of root powers
        subseries = [term]
        index = 0
        while index < len(series):
            if sorted(series[index]) == perm_id:
                subseries.append(series.pop(index))

            else:
                index += 1

        # Factor out constant multiple (all terms MUST have the same multiple)
        multiplicity = subseries.count(term)
        minimal_subseries = [*{*subseries}] # Remove duplicates
        assert all(subseries.count(other) == multiplicity for other in minimal_subseries)
        accumulate_series(reduced, distrib_const(multiplicity, reduce_prime_root_sum(minimal_subseries, esym, poly_deg)))

    return simplify_power_sum(reduced)

def tschirnhaus(poly_deg, trans_deg):
    trans_roots = []
    for root_id in range(poly_deg):
        # p + qx + rx^2 + sx^3 + ...
        trans_root = []
        for power in range(trans_deg + 1):
            trans_root.append((*tag_tuple(poly_deg, root_id, power), *tag_tuple(trans_deg, power)))

        trans_roots.append(trans_root)

    refactored = []
    for esym_index in range(poly_deg, 0, -1):
        sign = -1 if esym_index % 2 == 1 else 1
        groups = {}
        for combo in combos(trans_roots, esym_index):
            # Multiply out combination of transformed roots
            expansion = [*combo[0]]
            for other in combo[1:]:
                expansion = [
                    tuple(powa + powb for powa, powb in zip(left, right))
                    for left in expansion
                    for right in other
                ]

            # Sort by power-product tuples of transformation coefficients
            # i.e. {pq: ab+ac+ad+ae+bc+bd+be+cd+ce+de, p^2q: abc+abd+abe+..., ...}
            for term in expansion:
                root_prod, tcoeff_prod = term[:poly_deg], term[poly_deg:]
                if tcoeff_prod not in groups:
                    groups[tcoeff_prod] = []

                groups[tcoeff_prod].append(root_prod)

        # Reduce root power-product sums and distribute sign
        esym = [[*perms([1] * n + [0] * (poly_deg - n))] for n in range(1, poly_deg + 1)]
        coeff = []
        for tcoeff_prod, root_sum in groups.items():
            coeff.append([distrib_const(sign, reduce_root_sum(root_sum, esym, poly_deg)), tcoeff_prod])

        refactored.append(coeff)

    refactored.append(1)
    return refactored

#refactored = tschirnhaus(5, 4)
