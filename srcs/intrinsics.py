# 通过单应性矩阵计算内参
from srcs.utils import *


def generate_v_ij(H_stack, i, j):
    """Generate intrinsic orthogonality constraints. See Zhang pg. 6 for
       details.
    """
    N = H_stack.shape[0]
    v_ij = np.zeros((N, 6))
    v_ij[:, 0] = H_stack[:, 0, i] * H_stack[:, 0, j]
    v_ij[:, 1] = H_stack[:, 0, i] * H_stack[:, 1, j] + H_stack[:, 1, i] * H_stack[:, 0, j]
    v_ij[:, 2] = H_stack[:, 1, i] * H_stack[:, 1, j]
    v_ij[:, 3] = H_stack[:, 2, i] * H_stack[:, 0, j] + H_stack[:, 0, i] * H_stack[:, 2, j]
    v_ij[:, 4] = H_stack[:, 2, i] * H_stack[:, 1, j] + H_stack[:, 1, i] * H_stack[:, 2, j]
    v_ij[:, 5] = H_stack[:, 2, i] * H_stack[:, 2, j]
    return v_ij


def calc_intrinsics(Hs):
    N = len(Hs)
    H_stack = np.zeros((N, 3, 3))
    for idx, H in enumerate(Hs):
        H_stack[idx] = H
    # 注：实现时候我们-1
    v_00 = generate_v_ij(H_stack, i=0, j=0)
    v_01 = generate_v_ij(H_stack, i=0, j=1)
    v_11 = generate_v_ij(H_stack, i=1, j=1)

    V = np.zeros((2 * N, 6))
    V[:N] = v_01
    V[N:] = v_00 - v_11
    b = svd_solve(V)

    # B11, B12, B22, B13, B23, B33
    B0, B1, B2, B3, B4, B5 = b

    w = B0 * B2 * B5 - B1 ** 2 * B5 - B0 * B4 ** 2 + 2. * B1 * B3 * B4 - B2 * B3 ** 2
    d = B0 * B2 - B1 ** 2

    u0 = (B4 * B1 - B2 * B3) / d
    v0 = (B1 * B3 - B4 * B0) / d
    alpha = np.sqrt(w / (d * B0))
    beta = np.sqrt(w * B0 / (d ** 2))
    gamma = -B1 * np.sqrt(w / (B0 * d ** 2))
    gamma = 0
    K = np.array([[alpha, gamma, u0],
                   [0.,    beta,  v0],
                   [0.,    0.,    1.]])
    return K
