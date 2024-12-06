import numpy as np
from numpy import ndarray


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

    n = len(A)
    # for i in range(n):
    #     first_non_zero = -1
    #     for j in range(i + 1):
    #         if A[i, j] != 0:
    #             first_non_zero = j
    #             break
    #
    #     if first_non_zero != -1:
    #         values.extend(A[i, first_non_zero:i + 1])
    #
    #     pointers.append(len(values))

    for i in range(n):
        non_zero_indices = np.flatnonzero(A[i, :i + 1])  # 获取非零元素索引
        if len(non_zero_indices) > 0:
            first_non_zero = non_zero_indices[0]  # 第一个非零元素列索引
            values.extend(A[i, first_non_zero:i + 1])  # 添加从第一个非零到对角线的值
        pointers.append(len(values))  # 当前 values 的长度作为指针

    return values, pointers


K = np.array([
    [1, 1, 6, 0, 0],
    [1, 4, 1, 0, 0],
    [6, 1, 7, 1, 10],
    [0, 0, 1, 3, 1],
    [0, 0, 10, 1, 1]
])

v, p = skyline_storage(K)
print(v)
print(p)
