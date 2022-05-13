import os
from srcs.Calibrator import *

# 超参数
data_dir         = "data"
output_dir       = "outputs"
ends             = "jpg"
save_file_opencv = os.path.join(output_dir, "opencv.txt")
save_file_zhang  = os.path.join(output_dir, "zhang.txt")

points_per_row   = 9
points_per_col   = 6
square_size      = 10
show             = False



def main():
    # 01 读取图片路径
    files = [os.path.join(data_dir, file) for file in os.listdir(data_dir) if file.endswith(ends)]
    files.sort(key=lambda x: int(os.path.basename(x).split(".")[0]))

    # 02 调用opencv进行标定
    calibrator_opencv = CalibratorOpenCV()
    [img_size, points_world_xyz_all, points_pixel_xy_all] = calibrator_opencv.detectCorners(
        files, points_per_row, points_per_col, square_size, show)
    calibrator_opencv.calib(save_file_opencv)

    # 03 调用自己实现的张正友标定法
    calibrator_zhang = CalibratorZhang()
    calibrator_zhang.calib(img_size, points_world_xyz_all, points_pixel_xy_all)

if __name__ == '__main__':
    main()