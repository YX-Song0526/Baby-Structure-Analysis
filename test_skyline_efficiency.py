import skyline as sky
import numpy as np
from frame3D import s
import time

K = s.cal_K()
F = s.QnM

start_time = time.time()
u = sky.solve(K, F)
skyline_time = time.time() - start_time

start_time = time.time()
u1 = np.linalg.solve(K, F)
numpy_time = time.time() - start_time

# --------------------------- 结果输出与对比 ---------------------------
print(f"Skyline LDL 分解求解时间: {skyline_time:.6f}秒")
print(f"NumPy 内置求解时间: {numpy_time:.6f}秒")

