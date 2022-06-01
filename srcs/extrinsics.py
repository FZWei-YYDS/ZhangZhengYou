import numpy as np
from srcs.utils import re_orthogonalize


def recover_extrinsics(H, K):
    K_inv = np.linalg.inv(K)
    h0, h1, h2 = H[:, 0], H[:, 1], H[:, 2]
    lam = 1. / np.linalg.norm(np.dot(K_inv, h0))
    r0 = lam * np.dot(K_inv, h0)
    r1 = lam * np.dot(K_inv, h1)
    # 由于我们前面进行了归一化，因此直接进行向量乘法即可
    r2 = np.cross(r0, r1)

    t = lam * np.dot(K_inv, h2)
    # np.vstack 一行一行堆叠
    R = np.vstack((r0, r1, r2)).T
    # 由于数值计算的关系，R并非完全是正交矩阵，因此进行优化
    R = re_orthogonalize(R)
    E = np.hstack((R, t[:, np.newaxis]))
    return E