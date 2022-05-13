import os
import cv2
import numpy as np
from tqdm import tqdm

from srcs.homography import calc_homography

class CalibratorOpenCV(object):
    def detectCorners(self, files, points_per_row, points_per_col, square_size, show):
        self.points_world_xyz_all = []
        self.points_pixel_xy_all = []
        self.img_size = None

        # 建立角点的世界坐标
        points_world_xyz = np.zeros((points_per_row * points_per_col, 3), np.float32)
        points_world_xyz[:, :2] = np.mgrid[0: points_per_row, 0: points_per_col].T.reshape(-1, 2)
        points_world_xyz = points_world_xyz * square_size

        # 检点检测标准
        criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)
        for file in tqdm(files, "检测角点"):
            img = cv2.imread(file)
            if img is None:
                raise FileNotFoundError(file, "没有发现!")
            # 裁剪灰色图片
            if len(img.shape) == 2:
                gray = img
            else:
                gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
            # 检查图片尺寸是否一致
            if self.img_size is None:
                self.img_size = gray.shape
            else:
                assert gray.shape == self.img_size

            # 角点粗检测
            ret, points_pixel_xy = cv2.findChessboardCorners(gray, (points_per_row, points_per_col), None)
            # 精检测
            if ret:
                points_pixel_xy = cv2.cornerSubPix(gray, points_pixel_xy, (11, 11), (-1, -1), criteria)
                # 添加角点的像素坐标、世界坐标用于之后的标定
                self.points_world_xyz_all.append(points_world_xyz)
                self.points_pixel_xy_all.append(points_pixel_xy)
            else:
                print("未检测到角点:", file)

            if show:
                if len(img.shape) == 2:
                    # 转化为3通道的便于之后显示图像
                    img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
                img = cv2.drawChessboardCorners(img, (points_per_row, points_per_col), points_pixel_xy, ret)
                title = os.path.basename(file)
                cv2.imshow(title, img)
                cv2.moveWindow(title, 200, 100)
                cv2.waitKey(1000)
                cv2.destroyWindow(title)
        return [self.img_size, self.points_world_xyz_all, self.points_pixel_xy_all]

    def calib(self, save_file):
        self.error, self.mtx, self.dist, self.rvecs, self.tvecs = cv2.calibrateCamera(
            self.points_world_xyz_all,  # 世界坐标
            self.points_pixel_xy_all,  # 像素坐标
            self.img_size,  # 图像尺寸
            None, None
        )
        # 效果好坏评价
        with open(save_file, "w") as f:
            print("重投影误差:\n", self.error, "\n", file=f)
            print("相机内参:\n", self.mtx, "\n", file=f)
            print("相机畸变(k1,k2,p1,p2,k3):\n", self.dist, "\n", file=f)
            # 每个标定板（左上角角点）相对于相机原点的平移、旋转矩阵
            print("旋转矩阵:\n", self.rvecs, "\n", file=f)
            print("平移矩阵:\n", self.tvecs, "\n", file=f)

    def rectify(self, img):
        return cv2.undistort(img, self.mtx, self.dist)


class CalibratorZhang(CalibratorOpenCV):
    def calib(self, img_size, points_world_xyz_all, points_pixel_xy_all):
        # 01 计算初值
        # 1.1 计算单应性矩阵
        Hs = []
        for points_world_xyz, points_pixel_xy in zip(points_world_xyz_all, points_pixel_xy_all):
            points_world_xyz = points_world_xyz[:, :2]  # z=0默认
            points_pixel_xy = np.squeeze(points_pixel_xy)  # 删除多余轴
            # 计算单应性矩阵H
            H = calc_homography(points_world_xyz, points_pixel_xy)

            # 对单应性矩阵
            print(1)




