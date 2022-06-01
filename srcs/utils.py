import numpy as np


# 非齐次坐标 -> 齐次坐标
def to_homogeneous(x):
    return np.hstack((np.atleast_2d(x), np.ones((len(x), 1))))


# 齐次坐标 -> 非齐次坐标
def to_inhomogeneous(x):
    x = np.atleast_2d(x)  #
    x /= x[:, -1][:, np.newaxis]
    return x[:, :-1]


def svd_solve(A):
    """利用SVD求解齐次方程"""
    U, S, V_t = np.linalg.svd(A)
    # S已经从大到小排列了，方程解即奇异值最小的那列右奇异向量
    return V_t[np.argmin(S)]


def svd_solve2(A, b):
    """利用SVD求解非齐次方程"""
    U, S, V_t = np.linalg.svd(A)
    n = A.shape[-1]
    UnT = U[:, :n].T
    V = V_t.T
    r = np.linalg.matrix_rank(A)
    S_inv = np.linalg.inv(np.diag(S[:r]))
    x = np.matmul(np.matmul(np.matmul(V, S_inv),  UnT), b)
    return x


def pinv_solve(A, b):
    A_inv = np.linalg.pinv(A)
    return np.matmul(A_inv, b)


def re_orthogonalize(R):
    U, S, V_t = np.linalg.svd(R)
    return np.dot(U, V_t)


def save_result(save_file, error, A, dist):
    with open(save_file, "w") as f:
        print("保存标定结果到:", save_file)
        print("重投影误差:\n", error, "\n", file=f)
        print("重投影误差:\n", error, "\n")
        print("相机内参:\n", A, "\n", file=f)
        print("相机内参:\n", A, "\n")
        print("相机畸变(k1,k2,p1,p2,k3):\n", dist, "\n", file=f)
        print("相机畸变(k1,k2,p1,p2,k3):\n", dist, "\n")


if __name__ == '__main__':
    A = [[1, 2, 3]]
    A_hom = to_homogeneous(A)
    A_inhom = to_inhomogeneous(A_hom)
    print("齐次坐标：", A_hom)
    print("非齐次坐标：", A_inhom)
