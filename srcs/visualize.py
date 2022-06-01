import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from srcs.utils import to_homogeneous


def make_axis_publishable(ax, major_x, major_y, major_z):
    [t.set_va('center') for t in ax.get_yticklabels()]
    [t.set_ha('left') for t in ax.get_yticklabels()]
    [t.set_va('center') for t in ax.get_xticklabels()]
    [t.set_ha('right') for t in ax.get_xticklabels()]
    [t.set_va('center') for t in ax.get_zticklabels()]
    [t.set_ha('left') for t in ax.get_zticklabels()]

    ax.grid(False)
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.zaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    ax.xaxis._axinfo['tick']['inward_factor'] = 0
    ax.xaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.yaxis._axinfo['tick']['inward_factor'] = 0
    ax.yaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.zaxis._axinfo['tick']['inward_factor'] = 0
    ax.zaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.zaxis._axinfo['tick']['outward_factor'] = 0.4

    ax.xaxis.set_major_locator(MultipleLocator(major_x))
    ax.yaxis.set_major_locator(MultipleLocator(major_y))
    ax.zaxis.set_major_locator(MultipleLocator(major_z))


def camera_center(XYZ_all, extrinsic_matrices):
    XYZ_hom = to_homogeneous(XYZ_all[0])
    fig = plt.figure()
    ax = fig.add_subplot("111", projection="3d")

    make_axis_publishable(ax, 20, 20, 20)
    ax.set_title('Camera Center')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.set_xlim(-50, 50)
    ax.set_ylim(-50, 50)
    ax.set_zlim(0, 150)

    # From StackOverflow: https://stackoverflow.com/questions/39408794/python-3d-pyramid
    v = np.array([[-0.5, -0.5, 1],
                  [0.5, -0.5, 1],
                  [0.5, 0.5, 1],
                  [-0.5, 0.5, 1],
                  [0, 0, 0]])
    v = v * 8
    verts = [[v[0], v[1], v[4]], [v[0], v[3], v[4]],
             [v[2], v[1], v[4]], [v[2], v[3], v[4]],
             [v[0], v[1], v[2], v[3]]]

    ax.add_collection3d(Poly3DCollection(
        verts, facecolors='cyan', linewidths=0.5, edgecolors='r', alpha=.25))
    for E in extrinsic_matrices:
        XYZ_ext = np.dot(XYZ_hom, E.T)
        xs = XYZ_ext[:, 0]
        ys = XYZ_ext[:, 1]
        zs = XYZ_ext[:, 2]
        # 将每个角点绘制成面
        ax.plot_trisurf(xs, ys, zs)
    ax.invert_xaxis()
    plt.show()


def world_center(XYZ_all, extrinsic_matrices):
    XYZ_hom = to_homogeneous(XYZ_all[0])
    fig = plt.figure()

    ax = fig.add_subplot('111', projection='3d')

    make_axis_publishable(ax, 20, 20, 20)

    ax.set_title('World-Centric Extrinsics')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.set_xlim(-100, 100)
    ax.set_ylim(-100, 100)
    ax.set_zlim(-100, 20)

    # From StackOverflow: https://stackoverflow.com/questions/39408794/python-3d-pyramid
    v = np.array([[-0.5, -0.5, 1], [0.5, -0.5, 1], [0.5, 0.5, 1], [-0.5, 0.5, 1], [0, 0, 0]])

    v = v * 7
    x = 2
    xs = XYZ_hom[:, 0] * x
    ys = XYZ_hom[:, 1] * x
    zs = XYZ_hom[:, 2]
    ax.plot_trisurf(xs, ys, zs)
    v = to_homogeneous(v)

    for E in extrinsic_matrices:
        E = np.vstack((E, np.array([0., 0., 0., 1.])))
        E_inv = np.linalg.inv(E)
        E_inv = E_inv[:3]
        v_new = np.dot(v, E_inv.T)
        verts = [[v_new[0], v_new[1], v_new[4]], [v_new[0], v_new[3], v_new[4]],
                 [v_new[2], v_new[1], v_new[4]], [v_new[2], v_new[3], v_new[4]],
                 [v_new[0], v_new[1], v_new[2], v_new[3]]]
        ax.add_collection3d(Poly3DCollection(verts, facecolors='cyan', linewidths=0.5, edgecolors='r', alpha=.25))

    ax.invert_xaxis()
    plt.show()
