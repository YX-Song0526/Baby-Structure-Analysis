from space.systems import Frame3D
import numpy as np

s = Frame3D()

# 添加节点
s.add_node(0, 0, 0)
s.add_node(0, 2, 0)
s.add_node(2, 2, 0)
s.add_node(2, 0, 0)
s.add_node(0, 0, 3)
s.add_node(0, 2, 3)
s.add_node(2, 2, 3)
s.add_node(2, 0, 3)

# 添加单元
s.add_beam_column(1, 5, 69e9, 0.02, 2, 18, [0.05, 0.08])
s.add_beam_column(2, 6, 69e9, 0.02, 23, 18, [0.05, 0.08])
s.add_beam_column(3, 7, 69e9, 0.02, 2, 18, [0.05, 0.08])
s.add_beam_column(4, 8, 69e9, 0.02, 2, 18, [0.05, 0.08])
s.add_beam_column(8, 5, 69e9, 0.02, 2, 18, [0.05, 0.08])
s.add_beam_column(5, 6, 69e9, 0.02, 23, 18, [0.05, 0.08])
s.add_beam_column(6, 7, 69e9, 0.02, 2, 18, [0.05, 0.08])
s.add_beam_column(7, 8, 69e9, 0.02, 23, 18, [0.05, 0.08])

# 添加约束
s.add_fixed_sup(1, 2, 3, 4)

# 添加载荷
s.add_single_moment(6, Mx=0, My=2000000, Mz=-1500000)
# s.add_single_force(8, Fy=5000000)

# 求解
u = s.solve_una()
# q = s.solve_qnm()

# 打印结果
print('节点位移:',u)

# 可视化
s.plot_system(scale_max=1000)
