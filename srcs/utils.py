import numpy as np

# 非齐次坐标 -> 齐次坐标
def to_homogeneous(A):
    return np.hstack((np.atleast_2d(A), np.ones((len(A), 1))))

# 齐次坐标 -> 非齐次坐标
def to_inhomogeneous(A):
    A = np.atleast_2d(A)  #
    A /= A[:, -1][:, np.newaxis]
    return A[:, :-1]

def svd_solve(A):
    """利用SVD求解最小二乘法"""
    U, S, V_t = np.linalg.svd(A)
    idx = np.argmin(S)
    return V_t[idx]


if __name__ == '__main__':
    A = [[1, 2, 3]]
    A_hom = to_homogeneous(A)
    A_inhom = to_inhomogeneous(A_hom)
    print("齐次坐标：", A_hom)
    print("非齐次坐标：", A_inhom)