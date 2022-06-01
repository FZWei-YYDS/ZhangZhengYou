import numpy as np
from srcs import utils
from scipy.optimize import curve_fit


# 计算归一化的矩阵
def calc_N_matrix(A):
    if A.ndim != 2 or A.shape[-1] != 2:
        raise ValueError("必须为2D数据集")
    x, y = A[:, 0], A[:, 1]
    x_mean, y_mean = x.mean(), y.mean()
    x_var, y_var = x.var(), y.var()
    s_x, s_y = np.sqrt(2. / x_var), np.sqrt(2. / y_var)
    N_matrix = np.array([[s_x, 0., -s_x * x_mean],
                            [0., s_y, -s_y * y_mean],
                            [0,   0.,            1.]])
    return N_matrix  # 返回归一化矩阵


# 计算棋盘格的单应性矩阵初值
def calc_homography(XYs, uvs):
    N = XYs.shape[0]
    # 01 归一化数据，方便结果收敛
    Nxy = calc_N_matrix(XYs)
    Nuv = calc_N_matrix(uvs)

    # 02 将非齐次坐标转化为齐次坐标
    XY_hom = utils.to_homogeneous(XYs)
    uv_hom = utils.to_homogeneous(uvs)

    # 03 进行归一化
    # x' = N · x，但是代码中输入的是 x^T
    # 因此：矩阵乘法 (x')^T = (N · x)^T =  x^T N^
    XY_norm_hom = np.dot(XY_hom, Nxy.T)
    uv_norm_hom = np.dot(uv_hom, Nuv.T)

    # 04 解方程
    X, Y = XY_norm_hom[:, 0], XY_norm_hom[:, 1]
    u, v = uv_norm_hom[:, 0], uv_norm_hom[:, 1]

    M1 = np.zeros((N, 9))
    M1[:, 0] = -X; M1[:, 1] = -Y; M1[:, 2] = -1.
    M1[:, 6] = u * X; M1[:, 7] = u * Y; M1[:, 8] = u

    M2 = np.zeros((N, 9))
    M2[:, 3] = -X; M2[:, 4] = -Y; M2[:, 5] = -1.
    M2[:, 6] = v * X; M2[:, 7] = v * Y; M2[:, 8] = v

    M = np.zeros((N * 2, 9))
    M[: N], M[N: ] = M1, M2

    # 求解方程组，并得到
    h_norm = utils.svd_solve(M).reshape((3, 3))
    # 反归一化 $H = N_U^{ - 1} \cdot H' \cdot {N_X}$
    H = np.dot(np.dot(np.linalg.inv(Nuv), h_norm), Nxy)
    return H


def value_function_H(x_data, *params):
    """价值函数"""
    # 传入的P0参数的每个名字
    h11, h12, h13, h21, h22, h23, h31, h32, h33 = params
    N = x_data.shape[0] // 2
    # 世界坐标系
    X, Y = x_data[:N], x_data[N:]
    s = h31 * X + h32 * Y + h33
    u_proj = (h11 * X + h12 * Y + h13) / s
    v_proj = (h21 * X + h22 * Y + h23) / s
    uv_proj = np.zeros_like(x_data)
    # 为了之后求导方便
    uv_proj[:N], uv_proj[N:] = u_proj, v_proj
    return uv_proj


def jac_H(x_data, *params):
    h11, h12, h13, h21, h22, h23, h31, h32, h33 = params
    N = x_data.shape[0] // 2
    X, Y = x_data[:N], x_data[N:]

    # h11, h12, h13都是标量，所以直接用乘法
    s_u = h11 * X + h12 * Y + h13
    s_v = h21 * X + h22 * Y + h23
    w   = h31 * X + h32 * Y + h33
    w_2 = w ** 2

    J = np.zeros((N * 2, 9))
    J_u, J_v = J[:N, :], J[N:, :]
    J_u[:, 0] = X / w
    J_u[:, 1] = Y / w
    J_u[:, 2] = 1. / w
    J_u[:, 6] = (-s_u * X) / w_2  # 默认逐元素
    J_u[:, 7] = (-s_u * Y) / w_2
    J_u[:, 8] = -s_u / w_2

    J_v[:, 3] = X / w
    J_v[:, 4] = Y / w
    J_v[:, 5] = 1. / w
    J_v[:, 6] = (-s_v * X) / w_2
    J_v[:, 7] = (-s_v * Y) / w_2
    J_v[:, 8] = -s_v / w_2

    J[:N, :] = J_u
    J[N:, :] = J_v
    return J


# 执行非线性最小二乘法以细化线性单应性估计
def refine_homography(H, XYs, uvs):
    X, Y = XYs[:, 0], XYs[:, 1]
    u, v = uvs[:, 0], uvs[:, 1]
    N = X.shape[0]
    h0 = H.ravel()

    x_data = np.zeros(N * 2)
    x_data[:N], x_data[N:] = X, Y

    y_data = np.zeros(N * 2)
    y_data[:N], y_data[N:] = u, v

    # 01 非线性优化
    popt, pcov = curve_fit(f=value_function_H,  # 价值函数
                           xdata=x_data,  # X数据
                           ydata=y_data,  # Y数据
                           p0=h0,         # 初始值（参数组）
                           jac=jac_H,     # 价值函数的一阶导数（可选，不加scipy也会自动计算）
                           )
    H_refined = popt / popt[-1]
    H_refined = H_refined.reshape((3, 3))
    return H_refined
