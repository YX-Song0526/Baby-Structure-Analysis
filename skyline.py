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
               - values 是 Skyline 存储的非零值数组，
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
                break

        if first_nonzero != -1:
            values.extend(A[i, first_nonzero:i + 1])
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


def ldl_with_skyline(K: ndarray):
    KK = K.copy()
    K_sky, pointer, edge = skyline_storage(K)

    n = len(K)

    for j in range(1, n):
        # 计算g_ij
        for i in range(edge[j], j):
            sigma = 0
            for k in range(i):
                sigma += K_sky[pointer[i]-(i - k)] * K_sky[pointer[j]-(j - k)]
            K_sky[pointer[j] - (j - i)] -= sigma
            KK[j, i] -= sigma

        # 计算d_jj
        sigma = 0
        for k in range(edge[j], j):
            sigma += K_sky[pointer[j]-(j - k)] ** 2 / K_sky[pointer[k]]
        K_sky[pointer[j]] -= sigma
        KK[j, j] -= sigma

        # 计算L_ij
        for i in range(edge[j], j):
            K_sky[pointer[j] - (j - i)] /= K_sky[pointer[i]]
            KK[j, i] /= K_sky[pointer[i]]

    return K_sky, KK

# 示例矩阵
A = np.array([
    [4, 1, 2, 0],
    [1, 3, 0, 0],
    [2, 0, 6, 5],
    [0, 0, 5, 8]
])

# values, pointers, edge = skyline_storage(A)
# print("Values:", values)
# print("Pointers:", pointers)
# print("Edges:", edge)




