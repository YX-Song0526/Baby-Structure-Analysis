import numpy as np
from numpy import cos, sin


def little_mat(phi):
    """
    二维杆单元整体坐标下的单元刚度矩阵简化施工

    Args:
        phi: 偏转角（单位为弧度）

    Returns:
        4×4 矩阵

    """
    return np.array([[cos(phi) ** 2, cos(phi) * sin(phi), -cos(phi) ** 2, -cos(phi) * sin(phi)],
                     [cos(phi) * sin(phi), sin(phi) ** 2, -cos(phi) * sin(phi), -sin(phi) ** 2],
                     [-cos(phi) ** 2, -cos(phi) * sin(phi), cos(phi) ** 2, cos(phi) * sin(phi)],
                     [-cos(phi) * sin(phi), -sin(phi) ** 2, cos(phi) * sin(phi), sin(phi) ** 2]])


def tiny_mat(L):
    """
    一维EB梁单元单元刚度矩阵简化施工

    Args:
        L: 单元长度

    Returns:
        4×4 矩阵

    """
    return np.array([[12, 6 * L, -12, 6 * L],
                     [6 * L, 4 * L ** 2, -6 * L, 2 * L ** 2],
                     [-12, -6 * L, 12, -6 * L],
                     [6 * L, 2 * L ** 2, -6 * L, 4 * L ** 2]])


def big_mat(E, A, I, L):
    """
    二维钢架单元的局部单元刚度矩阵

    Args:
        E: 杨氏模量
        A: 截面积
        I: 惯性矩
        L: 单元长度

    Returns:
        6×6 单元刚度矩阵

    """
    return np.array([[E * A / L, 0, 0, -E * A / L, 0, 0],
                     [0, 12 * E * I / L ** 3, 6 * E * I / L ** 2, 0, -12 * E * I / L ** 3, 6 * E * I / L ** 2],
                     [0, 6 * E * I / L ** 2, 4 * E * I / L, 0, -6 * E * I / L ** 2, 2 * E * I / L],
                     [-E * A / L, 0, 0, E * A / L, 0, 0],
                     [0, -12 * E * I / L ** 3, -6 * E * I / L ** 2, 0, 12 * E * I / L ** 3, -6 * E * I / L ** 2],
                     [0, 6 * E * I / L ** 2, 2 * E * I / L, 0, -6 * E * I / L ** 2, 4 * E * I / L]])


def trans_mat_for_bar(phi):
    """
    二维干单元的坐标转换矩阵

    Args:
        phi: 偏转角（单位为弧度）

    Returns:
        4x4 的坐标转换矩阵

    """
    return np.array([[cos(phi), sin(phi), 0, 0],
                     [-sin(phi), cos(phi), 0, 0],
                     [0, 0, cos(phi), sin(phi)],
                     [0, 0, -sin(phi), cos(phi)]])


def trans_mat_for_frame(phi):
    """
    二维钢架单元坐标转换矩阵

    Args:
        phi: 偏转角（单位为弧度）

    Returns:
        6×6 的坐标转换矩阵

    """
    return np.array([[cos(phi), sin(phi), 0, 0, 0, 0],
                     [-sin(phi), cos(phi), 0, 0, 0, 0],
                     [0, 0, 1, 0, 0, 0],
                     [0, 0, 0, cos(phi), sin(phi), 0],
                     [0, 0, 0, -sin(phi), cos(phi), 0],
                     [0, 0, 0, 0, 0, 1]])
