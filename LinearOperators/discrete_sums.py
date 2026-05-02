import scipy
import numpy as np
from pyperclip import copy, paste

def sample_summ_op_diag(L, N):
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]
    I = np.eye(N)
    D = (I - np.concatenate((I[:, 1:], np.zeros((N, 1))), 1)) / dx
    S = np.linalg.inv(I - scipy.linalg.expm(-D))
    diag = np.flip(S - I, axis=1).diagonal()
    return diag / np.max(diag)

# f(n) = sum(g(k), 0 <= k <= n)
# f(n) - f(n - 1) = g(n); f(0) = 0
# (1 - e^-D)f = g; f(0) = 0
def solve_sum_approx(f, x0, x1, N):
    x = np.linspace(x0, x1, N)
    dx = x[1] - x[0]
    I = np.eye(N)
    D = (I - np.concatenate((I[:, 1:], np.zeros((N, 1))), 1)) / dx
    S = np.linalg.inv(I - scipy.linalg.expm(-D))
    f = f(x)
    F = np.dot(S, f)
    return locals()

def copy_to_desmos(data):
    copy("\\left[" + ", ".join(f"{d:.4f}" for d in data) + "\\right]")

if __name__ == "__main__":
    f = lambda x: np.exp(-x * x)
    env = solve_sum_approx(f, -15, 15, 2048)
    copy_to_desmos(env["F"])
