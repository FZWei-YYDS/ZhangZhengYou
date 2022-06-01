import cv2
from srcs.utils import *


# 参考：https://blog.csdn.net/qq_32998593/article/details/113063216
# matlab跟opencv的重投影计算方式有一些区别
def reproject_error(XYZ_all, uv_all, K, k, extrinsic_matrices):
    M = len(XYZ_all)
    N = XYZ_all[0].shape[0]
    total_error = 0.
    for idx, E in enumerate(extrinsic_matrices):
        uvs_obs = np.squeeze(uv_all[idx])
        uvs_obs = uvs_obs.astype(np.float64)
        uvs_pro = reproject_dist(K, k, E, XYZ_all[idx])
        for i in range(uvs_pro.shape[0]):
            uv_obs, uv_pro = uvs_obs[i, :], uvs_pro[i, :]
            error = cv2.norm(uv_obs, uv_pro, cv2.NORM_L2)
            total_error += error * error
    total_error = np.sqrt(total_error / (M * N))
    return total_error


def reproject(K, E, XYZs):
    XYZs_hom = to_homogeneous(XYZs)
    P = np.dot(K, E)  # 3 x 4
    uvs_hom = np.dot(XYZs_hom, P.T)
    return to_inhomogeneous(uvs_hom)


def reproject_dist(K, k, E, XYZs):
    XYZs_hom = to_homogeneous(XYZs)
    xys = to_inhomogeneous(np.dot(XYZs_hom, E.T))
    xys_dist = distort(k, xys)
    xys_dist_hom = to_homogeneous(xys_dist)
    uv_dist_hom = np.dot(xys_dist_hom, K.T)
    uv_dist = to_inhomogeneous(uv_dist_hom)
    return uv_dist


def distort(k, xys):
    xs, ys = xys[:, 0], xys[:, 1]
    rs = np.sqrt(xs ** 2 + ys ** 2)
    k0, k1 = k
    D = k0 * rs ** 2 + k1 * rs ** 4
    x_dist = xs * (1 + D)
    y_dist = ys * (1 + D)
    xys_dist = np.hstack((x_dist[:, np.newaxis], y_dist[:, np.newaxis]))
    return xys_dist
