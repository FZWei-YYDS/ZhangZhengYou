import os
from tqdm import tqdm
from srcs.homography import *
from srcs.intrinsics import *
from srcs.extrinsics import *
from srcs.distortion import *
from srcs.refinement import *
from srcs.model import *
from srcs.visualize import *


class CalibratorOpenCV(object):
    def detectCorners(self, files, points_per_row, points_per_col, square_size, show):
        """ 功能：检测棋盘格角点

        :param files: 需要标定的图片路径
        :param points_per_row: 每行的角点数量
        :param points_per_col: 每列的角点数量
        :param square_size:    棋盘格角点大小
        :param show:           是否显示角点

        :return:
        """
        XYZ_all = []
        uv_all = []
        img_size = None

        # 建立角点的世界坐标
        XYZs = np.zeros((points_per_row * points_per_col, 3), np.float32)
        XYZs[:, :2] = np.mgrid[0: points_per_row, 0: points_per_col].T.reshape(-1, 2)
        XYZs = XYZs * square_size

        # 检点检测标准
        criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)
        for file in tqdm(files, "检测角点"):
            img = cv2.imread(file)
            if img is None:
                raise FileNotFoundError(file, "没有发现文件!")
            # 裁剪灰色图片
            if len(img.shape) == 2:
                gray = img
            else:
                gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
            # 检查图片尺寸是否一致
            if img_size is None:
                img_size = gray.shape
            else:
                assert gray.shape == img_size
            # 角点粗检测（我们传入了points_per_row, points_per_col参数，因此可以自动区分棋盘格旋转）
            ret, uvs = cv2.findChessboardCorners(gray, (points_per_row, points_per_col), None)
            # 精检测
            if ret:
                uvs = cv2.cornerSubPix(gray, uvs, (11, 11), (-1, -1), criteria)
                # 添加角点的像素坐标、世界坐标用于之后的标定
                XYZ_all.append(XYZs)
                uv_all.append(uvs)
            else:
                print("未检测到角点:", file)

            if show:
                if len(img.shape) == 2:
                    img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR) # 转化为3通道的便于之后显示图像
                img = cv2.drawChessboardCorners(img, (points_per_row, points_per_col), uvs, ret)
                title = os.path.basename(file)
                cv2.imshow(title, img)
                cv2.moveWindow(title, 200, 100)
                cv2.waitKey(1000)
                cv2.destroyWindow(title)
        return [img_size, XYZ_all, uv_all]

    def calib(self, img_size, XYZ_all, uv_all):
        error, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(
            XYZ_all,   # 世界坐标
            uv_all,    # 像素坐标
            img_size,  # 图像尺寸
            None,      # 相机内参（初始）
            None,      # 相机畸变（初始）
            flags=cv2.CALIB_ZERO_TANGENT_DIST + cv2.CALIB_FIX_K3
            # CALIB_ZERO_TANGENT_DIST: 不估计切向畸变
            # CALIB_FIX_K3：不估计K3
            # 因为精确制作将导致很小的径向畸变。试图将参数拟合0会导致噪声干扰和数值不稳定
        )
        print("opencv的重投影误差:", error)
        return [error, mtx, dist, rvecs, tvecs]

    def rectify(self, img, mtx, dist):
        return cv2.undistort(img, mtx, dist)


class CalibratorZhang:
    def calib(self, XYZ_all, uv_all):
        # 01 计算初值
        # 1) 计算单应性矩阵
        Hs = []
        for XYZs, uvs in zip(XYZ_all, uv_all):
            # XYs：世界坐标系
            # uvs：像素坐标系
            XYs = XYZs[:, :2]      # z=0默认
            uvs = np.squeeze(uvs)  # 删除多余轴

            # 01 DLT计算单应性矩阵初值
            H = calc_homography(XYs, uvs)
            # 02 使用重投影误差非线性优化单应性矩阵
            H = refine_homography(H, XYs, uvs)
            Hs.append(H)
        # 2) 计算相机内参
        K = calc_intrinsics(Hs)

        # 3) 计算相机外参
        extrinsic_matrices = []
        for idx, H in enumerate(Hs):
            # E: [R:t] 3 x 4
            E = recover_extrinsics(H, K)
            extrinsic_matrices.append(E)

        # 4) 计算镜头畸变（相机内参、外参不变的情况下）
        k = calculate_lens_distortion(XYZ_all, uv_all, K, extrinsic_matrices)

        # 5) 初值投影误差
        error1 = reproject_error(XYZ_all, uv_all, K, k, extrinsic_matrices)
        print("重投影误差（线性）:", error1)

        # 6) 非线性优化所有参数
        K_opt, k_opt, extrinsic_matrices_opt = refine_all_parameters(XYZ_all, uv_all, K, k, extrinsic_matrices)

        error2 = reproject_error(XYZ_all, uv_all, K_opt, k_opt, extrinsic_matrices_opt)
        print("重投影误差（非线性）:", error2)

        # 7) 查看标定结果
        camera_center(XYZ_all, extrinsic_matrices_opt)
        world_center(XYZ_all, extrinsic_matrices_opt)

        return error2, K_opt, k_opt, extrinsic_matrices_opt