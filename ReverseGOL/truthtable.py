from more_itertools import distinct_permutations
import numpy as np

def karnaugh_reduce(table):
    kmap = []
    uncovered = np.ones(table.shape, dtype=np.bool)
    table_view, uncovered_view = table.view(), uncovered.view()
    for index, value in np.ndenumerate(table):
        if not (value and uncovered[index]):
            continue

        for dims in range(len(table.shape), 0, -1):
            best_slice = slice(None)
            max_new_covered = 0
            for axis_vars in distinct_permutations([i < dims for i in range(len(table.shape))]):
                block_slice = tuple(slice(None) if axis_var else subindex for subindex, axis_var in zip(index, axis_vars))
                if np.all(table_view[block_slice]):
                    new_covered = np.sum(uncovered_view[block_slice])
                    if new_covered > max_new_covered:
                        best_slice = block_slice
                        max_new_covered = new_covered

            if max_new_covered > 0:
                kmap.append(best_slice)
                uncovered[best_slice] = np.False_
                break

    return kmap

def karnaugh_reduce_ex(table, dont_cares):
    kmap = []
    uncovered = np.ones(table.shape, dtype=np.bool)
    max_table = np.logical_or(table, dont_cares)
    max_table_view, uncovered_view = max_table.view(), uncovered.view()
    for index, value in np.ndenumerate(table):
        if not (value and uncovered[index]):
            continue

        for dims in range(len(table.shape), 0, -1):
            best_slice = slice(None)
            max_new_covered = 0
            for axis_vars in distinct_permutations([i < dims for i in range(len(table.shape))]):
                block_slice = tuple(slice(None) if axis_var else subindex for subindex, axis_var in zip(index, axis_vars))
                if np.all(max_table_view[block_slice]):
                    new_covered = np.sum(uncovered_view[block_slice])
                    if new_covered > max_new_covered:
                        best_slice = block_slice
                        max_new_covered = new_covered

            if max_new_covered > 0:
                kmap.append(best_slice)
                uncovered[best_slice] = np.False_
                break

    return kmap
