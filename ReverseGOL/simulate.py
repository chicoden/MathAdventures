import time

rule = [
    [0, 0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 0, 0, 0, 0, 0]
]

def print_grid(grid):
    for row in grid:
        print(*("x" if cell else " " for cell in row))

def gol_step_forward(grid):
    next_grid = [[0 for cell in row] for row in grid]
    for y in range(len(grid)):
        for x in range(len(grid[y])):
            live_neighbors = 0
            for neighbor_y in range(max(0, y - 1), min(y + 2, len(grid))):
                for neighbor_x in range(max(0, x - 1), min(x + 2, len(grid[y]))):
                    if neighbor_x != x or neighbor_y != y:
                        live_neighbors += grid[neighbor_y][neighbor_x]

            next_grid[y][x] = rule[grid[y][x]][live_neighbors]

    return next_grid

def gol_simulate(grid, time_step=1, steps=-1):
    i = 0
    print_grid(grid)
    while i != steps:
        print("-" * len(grid[-1]) * 2)
        time.sleep(time_step)
        print_grid(grid := gol_step_forward(grid))
        i += 1

if __name__ == "__main__":
    gol_simulate([
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1, 1, 1],
        [0, 1, 1, 0, 0, 0, 1, 0],
        [0, 1, 0, 1, 0, 0, 0, 0],
        [0, 0, 1, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ])
