from truthtable import karnaugh_reduce
import numpy as np

rule = [
    [0, 0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 0, 0, 0, 0, 0]
]

grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0],
    [0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0],
    [0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0],
    [0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0],
    [0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

grid_width = len(grid[0])
grid_height = len(grid)

def kmap_to_maxterms(kmap):
    return [[(index, axis_slice) for index, axis_slice in enumerate(block_slice) if axis_slice != slice(None)] for block_slice in kmap]

def simplify_maxterms(maxterms, fixed_vars):
    new_maxterms = []
    for maxterm in maxterms:
        if any(fixed_vars[index] != negate for index, negate in maxterm if fixed_vars[index] != -1):
            continue

        new_maxterm = [(index, negate) for index, negate in maxterm if fixed_vars[index] == -1]
        if len(new_maxterm) > 0:
            new_maxterms.append(new_maxterm)

    return new_maxterms

def simplify_complement_maxterms(inv_maxterms, maxterms, fixed_vars):
    return (simplify_maxterms(inv_maxterms, fixed_vars), simplify_maxterms(maxterms, fixed_vars))

print("Generating update truthtable...")
truthtable = np.zeros([2] * 9, dtype=np.bool)
for index, value in np.ndenumerate(truthtable):
    truthtable[index] = rule[index[4]][sum(index) - index[4]]
print("Done.")

print("Computing Karnaugh maps...")
alive_kmap = karnaugh_reduce(~truthtable)
dead_kmap = karnaugh_reduce(truthtable)
print("Done.")

print("Creating expression templates...")
central_maxterms = (kmap_to_maxterms(dead_kmap), kmap_to_maxterms(alive_kmap))
corner_tl_maxterms = simplify_complement_maxterms(*central_maxterms, (0, 0, 0, 0, -1, -1, 0, -1, -1))
corner_tr_maxterms = simplify_complement_maxterms(*central_maxterms, (0, 0, 0, -1, -1, 0, -1, -1, 0))
corner_bl_maxterms = simplify_complement_maxterms(*central_maxterms, (0, -1, -1, 0, -1, -1, 0, 0, 0))
corner_br_maxterms = simplify_complement_maxterms(*central_maxterms, (-1, -1, 0, -1, -1, 0, 0, 0, 0))
edge_t_maxterms = simplify_complement_maxterms(*central_maxterms, (0, 0, 0, -1, -1, -1, -1, -1, -1))
edge_b_maxterms = simplify_complement_maxterms(*central_maxterms, (-1, -1, -1, -1, -1, -1, 0, 0, 0))
edge_l_maxterms = simplify_complement_maxterms(*central_maxterms, (0, -1, -1, 0, -1, -1, 0, -1, -1))
edge_r_maxterms = simplify_complement_maxterms(*central_maxterms, (-1, -1, 0, -1, -1, 0, -1, -1, 0))

neighborhood_lookup = {
    0b111_111_111: central_maxterms,
    0b110_110_000: corner_tl_maxterms,
    0b011_011_000: corner_tr_maxterms,
    0b000_110_110: corner_bl_maxterms,
    0b000_011_011: corner_br_maxterms,
    0b111_111_000: edge_t_maxterms,
    0b000_111_111: edge_b_maxterms,
    0b110_110_110: edge_l_maxterms,
    0b011_011_011: edge_r_maxterms,
}
print("Done.")

print("Enumerating conditions...")
conditions = []
for y, row in enumerate(grid):
    for x, cell in enumerate(row):
        neighbors = [-1, -1, -1, -1, -1, -1, -1, -1, -1]
        neighborhood = 0b000_000_000
        for other_y in range(max(0, y - 1), min(grid_height, y + 2)):
            for other_x in range(max(0, x - 1), min(grid_width, x + 2)):
                neighborhood_index = (other_y - y + 1) * 3 + (other_x - x + 1)
                neighbors[neighborhood_index] = other_y * grid_width + other_x + 1
                neighborhood |= 1 << neighborhood_index

        maxterms = neighborhood_lookup[neighborhood][cell]
        conditions.extend([-neighbors[index] if negate else neighbors[index] for index, negate in maxterm] for maxterm in maxterms)
print("Done.")

print("Writing CNF file...")
with open("reverse_gol.cnf", "w") as file:
    file.write("p cnf {0} {1}\n".format(grid_width * grid_height, len(conditions)))
    for condition in conditions:
        file.write(" ".join(map(str, condition)) + " 0\n")
print("Done.")
