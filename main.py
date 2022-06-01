from srcs.calibrator import *
from srcs.utils import save_result


# 超参数
data_dir         = "data"
output_dir       = "outputs"
ends             = "jpg"
save_file_opencv = os.path.join(output_dir, "opencv.txt")
save_file_zhang  = os.path.join(output_dir, "zhang.txt")

points_per_row   = 9         # 每行棋盘格角点数量
points_per_col   = 6         # 每列棋盘格角点数量
square_size      = 10.       # 单个棋盘格的物理尺寸大小(mm)
show             = False      # 是否展示棋盘格检测角点


def main():
    # 01 读取标定图片路径
    files = [os.path.join(data_dir, file) for file in os.listdir(data_dir) if file.endswith(ends)]
    files.sort(key=lambda x: int(os.path.basename(x).split(".")[0]))

    # 02 调用opencv进行标定
    calibrator_opencv = CalibratorOpenCV()
    [img_size, XYZ_all, uv_all] = calibrator_opencv.detectCorners(
        files, points_per_row, points_per_col, square_size, show)

    error_opencv, A_opencv, k_opencv, rvecs, tvecs = calibrator_opencv.calib(
        img_size, XYZ_all, uv_all)

    # 03 调用自己实现的张正友标定法
    calibrator_zhang = CalibratorZhang()
    error_zhang, A_zhang, k_zhang, rts = calibrator_zhang.calib(XYZ_all, uv_all)

    # 04 移除图像畸变
    img = cv2.imread(files[0], 0)
    h, w = img.shape[:2]
    # 畸变校正，提前计算映射矩阵
    map_opencv_x, map_opencv_y = cv2.initUndistortRectifyMap(
        A_opencv, k_opencv, None, A_opencv, (w, h), 5)

    k_zhang = np.array([k_zhang[0], k_zhang[1], 0., 0., 0.])
    map_zhang_x, map_zhang_y = cv2.initUndistortRectifyMap(
         A_zhang, k_zhang, None, A_zhang, (w, h), 5)

    # 校正后的图片相当于没有畸变了，直接用相机内参即可，无需k1, k2
    img_opencv = cv2.remap(img, map_opencv_x, map_opencv_y, cv2.INTER_LINEAR)
    img_zhang = cv2.remap(img, map_zhang_x, map_zhang_y, cv2.INTER_LINEAR)

    # 05 展示结果
    cv2.imshow("img", img)
    cv2.imshow("img_opencv", img_opencv)
    cv2.imshow("img_zhang", img_zhang)
    cv2.imshow("dif:img - img_zhang", img - img_zhang)
    cv2.imshow("dif:img_opencv - img_zhang", img_opencv - img_zhang)

    # 06 写入标定结果
    save_result(save_file_opencv, error_opencv, A_opencv, k_opencv)
    save_result(save_file_zhang, error_zhang, A_zhang, k_zhang)

    cv2.waitKey(0)
    cv2.destroyAllWindows()


if __name__ == '__main__':
    main()
