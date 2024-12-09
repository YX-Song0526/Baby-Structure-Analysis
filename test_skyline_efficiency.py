import solver.skyline as sky
import numpy as np
from frame3D import s
import time

print(len(s.cal_K_total()))
# 获取输入矩阵K和载荷向量F
K = s.cal_K_total()[:40, :40]  # 计算刚度矩阵K
F = s.FnM[:40]  # 载荷向量

# --------------------------- 重复多次求解 ---------------------------
num_iterations = 100  # 重复求解次数
skyline_times = []
numpy_times = []

for _ in range(num_iterations):
    # --------------------------- 使用 Skyline 算法求解 ---------------------------
    start_time = time.time()
    u = sky.solve(K, F)  # 使用Skyline LDL分解求解
    skyline_times.append((time.time() - start_time) * 1000)  # 转换为毫秒

    # --------------------------- 使用 NumPy 的内置求解器 ---------------------------
    start_time = time.time()
    u1 = np.linalg.solve(K, F)  # 使用 NumPy 的内置函数求解线性方程
    numpy_times.append((time.time() - start_time) * 1000)  # 转换为毫秒

# 计算平均时间
avg_skyline_time = np.mean(skyline_times)
avg_numpy_time = np.mean(numpy_times)

# --------------------------- 结果输出与对比 ---------------------------
print(f"Skyline LDL 分解平均求解时间: {avg_skyline_time:.6f} 毫秒")
print(f"NumPy 内置求解平均时间: {avg_numpy_time:.6f} 毫秒")

# 输出两个解的差异，检查是否接近（如果两个解非常接近，则说明算法基本正确）
error = np.linalg.norm(u - u1)  # 计算最后一次解的L2范数误差
print(f"解的误差 (L2范数): {error:.6e}")

# 可选：判断解是否基本相同，通常会设置一个容差
tolerance = 1e-6
if error < tolerance:
    print("两个解非常接近，算法正确！")
else:
    print("两个解差异较大，请检查代码。")
