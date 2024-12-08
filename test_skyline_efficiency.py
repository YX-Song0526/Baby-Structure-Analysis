import numpy as np
import time
from test_skyline_ldl import ldl, ldl_modified, ldl_modified_with_skyline  # 导入你在 ldl.py 中定义的函数
from truss2D import s1
from frame3D import s


def test_ldl_performance(A):
    # 你提供的矩阵 A 将用于性能测试
    print(f"Testing matrix of size: {A.shape[0]}x{A.shape[1]}")

    # 测试 ldl
    start_time = time.time()
    l, d = ldl(A)
    ldl_time = time.time() - start_time

    # 测试 ldl_modified
    start_time = time.time()
    l_modified, d_modified = ldl_modified(A)
    ldl_modified_time = time.time() - start_time

    # 测试 ldl_modified_with_skyline
    start_time = time.time()
    l_modified_skyline, d_modified_skyline = ldl_modified_with_skyline(A)
    ldl_modified_skyline_time = time.time() - start_time

    # 输出结果
    print(f"  LDL time: {ldl_time:.4f} seconds")
    print(f"  LDL_modified time: {ldl_modified_time:.4f} seconds")
    print(f"  LDL_modified_Skyline time: {ldl_modified_skyline_time:.4f} seconds")

    # 返回结果
    return {
        'size': A.shape[0],
        'ldl_time': ldl_time,
        'ldl_modified_time': ldl_modified_time,
        'ldl_modified_skyline_time': ldl_modified_skyline_time,
    }

# 调用性能测试
if __name__ == "__main__":
    # 指定一个矩阵进行测试
    A = s.cal_K()

    # 进行性能测试
    test_results = test_ldl_performance(A)
    print(test_results)
