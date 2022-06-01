from scipy.optimize import curve_fit
from srcs.model import *


def pack_params(K, k, extrinsic_matrices):
    packed_params = []
    alpha, beta, gamma, u0, v0 = K[0, 0], K[1, 1], K[0, 1], K[0, 2], K[1, 2]
    k0, k1 = k
    # 不顾及gamma
    # a = [alpha, beta, gamma, u0, v0, k0, k1]
    a = [alpha, beta, u0, v0, k0, k1]
    packed_params.extend(a)
    for E in extrinsic_matrices:
        R = E[:3, :3]
        t = E[:, 3]
        rodrigues = cv2.Rodrigues(R)[0]
        rho_x, rho_y, rho_z = rodrigues
        t_x, t_y, t_z = t
        w = [rho_x, rho_y, rho_z, t_x, t_y, t_z]
        packed_params.extend(w)
    return np.array(packed_params, dtype=np.float64)


def unpack_refinement_params(params):
    intrinsics = params[:6]
    alpha, beta, u0, v0, k0, k1 = intrinsics
    K = np.array([[alpha, 0, u0],
                  [0.,    beta,  v0],
                  [0.,    0.,    1.]])
    k = np.array([k0, k1])
    # 我们使用6个参数r1,r2,r3,t1,t2,t3
    extrinsic_matrices = []
    for i in range(6, len(params), 6):
        rho_x, rho_y, rho_z, t_x, t_y, t_z = params[i: i + 6]
        R = cv2.Rodrigues(np.array([rho_x, rho_y, rho_z]))[0]
        t = np.array([t_x, t_y, t_z])
        E = np.zeros((3, 4))
        E[:3, :3] = R
        E[:, 3] = t
        extrinsic_matrices.append(E)
    return K, k, extrinsic_matrices


def refine_all_parameters(XYZ_all, uv_all, K, k, extrinsic_matrices):
    M = len(XYZ_all)
    N = XYZ_all[0].shape[0]
    p = pack_params(K, k, extrinsic_matrices)
    XY = XYZ_all[0][:, :2]
    X, Y = XY[:, 0], XY[:, 1]
    # 世界坐标
    X_data = np.zeros(2 * M * N)
    for i in range(M):
        X_data[i * N: (i + 1) * N] = X
        X_data[M * N + i * N: M * N + (i + 1) * N] = Y

    # 像素坐标
    obs_u, obs_v = [], []
    for i, uvs in enumerate(uv_all):
        uvs = np.squeeze(uvs)
        us, vs = uvs[:, 0], uvs[:, 1]
        obs_u.append(us)
        obs_v.append(vs)
    obs_u = np.hstack(obs_u)
    obs_v = np.hstack(obs_v)
    Y_data = np.hstack((obs_u, obs_v))

    # 通过LM方法最小化重投影误差（这里没有雅克比矩阵，scipy默认通过数值的方法计算）
    popt, pconv = curve_fit(f_refine_all, X_data, Y_data, p)
    K_opt, k_opt, extrinsic_matrices_opt = unpack_refinement_params(popt)
    return K_opt, k_opt, extrinsic_matrices_opt



# 前向传播过程
def f_refine_all(X_data, *params):
    K, k, extrinsic_matrices = unpack_refinement_params(params)
    M = len(extrinsic_matrices)
    N = X_data.shape[0] // (2 * M)
    # 世界坐标系,只需要取前面N列（每个棋盘格都是一样的）
    X = X_data[:M * N][:N]
    Y = X_data[M * N:][:N]
    Z = np.zeros_like(X)
    XYZs = np.vstack((X, Y, Z)).T

    us_pro = np.zeros(M * N)
    vs_pro = np.zeros(M * N)
    for idx, E in enumerate(extrinsic_matrices):
        uv_pro = reproject_dist(K, k, E, XYZs)
        us_pro[idx * N: (idx + 1) * N] = uv_pro[:, 0]
        vs_pro[idx * N: (idx + 1) * N] = uv_pro[:, 1]
    uvs_pro = np.zeros(2 * M * N)
    uvs_pro[:M * N] = us_pro
    uvs_pro[M * N:] = vs_pro
    return uvs_pro
