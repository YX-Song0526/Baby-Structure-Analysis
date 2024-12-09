import numpy as np
from numpy import ndarray
from truss2D import s1
import pandas as pd


def skyline_storage(A: ndarray):
    """
    将矩阵 A 转换为 Skyline 存储格式。

    Args:
        A (np.ndarray): 输入的稀疏对称矩阵。

    Returns:
        tuple: (values, pointers)
               - values 是轮廓线内元素数组，
               - pointers 是对角线元素在 values 中的索引数组。
    """
    values = []
    pointers = []
    edge = []

    n = len(A)
    for i in range(n):
        first_nonzero = -1
        for j in range(i + 1):
            if A[i, j] != 0:
                first_nonzero = j
                # print(first_nonzero)
                break

        if first_nonzero != -1:
            values.extend(A[i, first_nonzero:i + 1])
            # print(A[i, first_nonzero:i + 1])
            edge.append(first_nonzero)

        pointers.append(len(values) - 1)

    return values, pointers, edge


def get_semi_bandwidth(A: ndarray):
    """获取半带宽数组"""
    semi_bandwidth = []
    for i in range(len(A)):
        for j in range(i + 1):
            if A[i, j] != 0:
                semi_bandwidth.append(i - j)
                break

    return semi_bandwidth


def skyline_to_full_matrix(values, pointers):
    """
    将 Skyline 压缩存储的矩阵转换为普通的矩阵

    Args:
        values: 轮廓线内元素数组
        pointers: 对角线元素在 values 中的索引数组

    Returns: n维矩阵 (np.ndarray)

    """



def ldl_modified_with_skyline(K):
    n = len(K)

    K_sky, pointer, edge = skyline_storage(K)

    for j in range(n):
        # 求g_ij
        for i in range(edge[j], j):
            for k in range(max(edge[i], edge[j]), i):
                K_sky[pointer[j] - (j - i)] -= K_sky[pointer[i] - (i - k)] * K_sky[pointer[j] - (j - k)]

        # 求D_ij
        for k in range(edge[j], j):
            K_sky[pointer[j]] -= K_sky[pointer[j] - (j - k)] ** 2 / K_sky[pointer[k]]

        # 求L_ij
        for i in range(edge[j], j):
            K_sky[pointer[j] - (j - i)] = K_sky[pointer[j] - (j - i)] / K_sky[pointer[i]]

    return K_sky
