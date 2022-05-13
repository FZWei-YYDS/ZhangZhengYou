import cv2
import numpy as np
from srcs import utils


# 计算归一化的矩阵
def calc_norm_matrix(A):
    if A.ndim != 2 or A.shape[-1] != 2:
        raise ValueError("必须为2D数据集")
    x, y = A[:, 0], A[:, 1]
    x_mean, y_mean = x.mean(), y.mean()
    x_var, y_var = x.var(), y.var()
    s_x, s_y = np.sqrt(2. / x_var), np.sqrt(2. / y_var)
    norm_matrix = np.array([[s_x, 0., -s_x * x_mean],
                            [0., s_y, -s_y * y_mean],
                            [0,   0.,            1.]])
    return norm_matrix


# 计算棋盘格的单应性矩阵初值
def calc_homography(points_world_xy, points_pixel_xy):
    N = points_world_xy.shape[0]
    # 归一化数据，方便结果收敛
    norm_matrix_world = calc_norm_matrix(points_world_xy)
    norm_matrix_pixel = calc_norm_matrix(points_pixel_xy)

    # 转化为齐次坐标
    world_hom = utils.to_homogeneous(points_world_xy)
    pixel_hom = utils.to_homogeneous(points_pixel_xy)

    # 进行归一化
    world_hom_norm = np.dot(world_hom, norm_matrix_world.T)
    pixel_hom_norm = np.dot(pixel_hom, norm_matrix_pixel.T)

    X, Y = world_hom_norm[:, 0], world_hom_norm[:, 1]
    u, v = pixel_hom_norm[:, 0], pixel_hom_norm[:, 1]

    # 齐次约束矩阵
    M = np.zeros((N * 2, 9))

    u_component = np.zeros((N, 9))
    u_component[:, 0] = -X
    u_component[:, 1] = -Y
    u_component[:, 2] = -1.
    u_component[:, 6] = u * X
    u_component[:, 7] = u * Y
    u_component[:, 8] = u

    v_component = np.zeros((N, 9))
    v_component[:, 3] = -X
    v_component[:, 4] = -Y
    v_component[:, 5] = -1.
    v_component[:, 6] = v * X
    v_component[:, 7] = v * Y
    v_component[:, 8] = v

    # 为了方便，把
    M[: N] = u_component
    M[N: ] = v_component

    h_norm = utils.svd_solve(M).reshape((3, 3))
    H = np.dot(np.dot(np.linalg.inv(norm_matrix_pixel), h_norm), norm_matrix_world)
    return H




# 执行非线性最小二乘法以细化线性单应性估计
def refine_homography(H, points_world_xy, points_pixel_xy):
    X, Y = points_world_xy[:, 0], points_world_xy[:, 1]
    u, v = points_pixel_xy[:, 0], points_pixel_xy[:, 1]

    N = X.shape[0]
    h0 = H.ravel()

    x_data = np.zeros(N * 2)
    x_data[: N] = X
    x_data[N: ] = Y

    y_data = np.zeros(N * 2)
    y_data[: N] = u
    y_data[N: ] = v

    # 使用Levenberg-Marquardt细化线性单应性估计


