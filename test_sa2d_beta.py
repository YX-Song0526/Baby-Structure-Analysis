from sa2d_beta.systems import Frame2D
from sa2d_beta.section import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

s = Frame2D()

# 创建材料
steel = Material(E=69e9)

# 创建截面
section_1 = Section(material=steel, shape=Circle(0.05))

# 添加节点
s.add_node(-4., 0.)  # node1
s.add_node(0., 0.)  # node2
s.add_node(2., 0.)  # node3
s.add_node(0., 2.)  # node4
s.add_node(0., 4.)  # node5

# 添加钢架并指定截面
s.add_element(1, 5, section_1)  # e1
s.add_element(2, 4, section_1)  # e2
s.add_element(4, 5, section_1)  # e3
s.add_element(4, 3, section_1)  # e4

# 添加边界条件
s.add_fixed_sup(1, 2, 3)

# 添加力
s.add_single_force(5, Fx=50000)

U = s.solve_disp()
max_stress = s.get_max_stress()
print(U)
print(max_stress)

# 计算最大应力
max_stress_list = []
R_range = np.arange(0.01, 0.1, 0.005)
for R in R_range:
    for element in s.elements:
        element.update_shape_params(R=R)
    max_stress = s.get_max_stress()
    max_stress = max_stress / 1e6
    max_stress_list.append(max_stress)


# 定义你的拟合函数，例如多项式函数
def polynomial(x, a, b):
    return a * (1 / x ** 4) + b * (1 / x ** 2)


# 使用curve_fit进行拟合
params, covariance = curve_fit(polynomial, R_range, max_stress_list)

# 拟合结果的参数
a, b = params
print(f"拟合参数: a={a}, b={b}")

# 绘制结果
x_fit = np.linspace(min(R_range), max(R_range), 100)
y_fit = polynomial(x_fit, *params)

plt.plot(x_fit, y_fit, label='拟合曲线', color='blue')
plt.scatter(R_range, max_stress_list, alpha=0.5)
plt.show()
