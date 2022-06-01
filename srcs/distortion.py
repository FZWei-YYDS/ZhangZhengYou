from srcs.utils import *
from srcs.model import reproject


def calculate_lens_distortion(XYZ_all, uv_all, K, extrinsic_matrices):
    # 多少副标定图片
    M = len(uv_all)
    N = uv_all[0].shape[0]

    u0, v0 = K[0, 2], K[1, 2]

    # 01 计算r_{i,j}
    r = np.zeros(2 * M * N)

    for idx, E in enumerate(extrinsic_matrices):
        # 将坐标转换为齐次坐标 54 x 4
        XYZ_hom = to_homogeneous(XYZ_all[idx])
        xy_norm = to_inhomogeneous(np.dot(XYZ_hom, E.T))
        x_norm, y_norm = xy_norm[:, 0], xy_norm[:, 1]
        r_i = np.sqrt(x_norm ** 2 + y_norm ** 2)
        r[idx * N: (idx + 1)  * N] = r_i
    # 我们对方程进行了重新排布
    r[M * N:] = r[:M * N]

    # 02 畸变后的像素点坐标  u_{ij}_dist, v_{ij}_dist
    uv_dists = np.zeros(2 * M * N)
    u_dists, v_dists = np.zeros(M * N), np.zeros(M * N)
    for idx, uvs in enumerate(uv_all):
        uvs = np.squeeze(uvs)
        u_dists[idx * N: (idx + 1) * N] = uvs[:, 0]
        v_dists[idx * N: (idx + 1) * N] = uvs[:, 1]
    uv_dists[: M * N] = u_dists
    uv_dists[M * N:] = v_dists

    # 03 理想无畸变像素坐标（根据相机参数计算的重投影像素坐标）
    uv_reps = np.zeros(2 * M * N)
    u_reps, v_reps = np.zeros(M * N), np.zeros(M * N)
    for idx, E in enumerate(extrinsic_matrices):
        uv_rep = reproject(K, E, XYZ_all[idx])
        u_reps[idx * N: (idx + 1) * N] = uv_rep[:, 0]
        v_reps[idx * N: (idx + 1) * N] = uv_rep[:, 1]
    uv_reps[:M * N] = u_reps
    uv_reps[M * N:] = v_reps

    # 04 构建方程组
    p = np.zeros(2 * M * N)
    p[:M * N] = u_reps - u0
    p[M * N:] = v_reps - v0

    D = np.zeros((2 * M * N, 2))
    D[:, 0] = p * r ** 2
    D[:, 1] = p * r ** 4

    d = uv_dists - uv_reps

    # 使用SVD求解方程
    ks1 = svd_solve2(D, d)
    # 使用伪逆求解方程
    # ks2 = pinv_solve(D, d)
    # 两者的差值很小
    # dif = ks1 - ks2
    return ks1
