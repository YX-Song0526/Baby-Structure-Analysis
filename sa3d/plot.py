import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from sa3d.elements import Node, Bar, BeamColumn
from typing import Union, List
import numpy as np


def plot_system(elements: List[Union[Bar, BeamColumn]],
                nodes: List[Node],
                U: np.ndarray,
                n_dof_per_node: int,
                initial_scale=1.0,
                max_scale=10000.0):
    """
    绘制结构系统的原始和变形状态。

    Args:
    - elements: 包含元素（如 Bar 或 BeamColumn）的列表。
    - nodes: 包含节点（Node 对象）的列表。
    - displacements: 节点位移向量。
    - n_dof_per_node: 每个节点的自由度数量（通常为3，表示x, y, z位移）。
    - initial_scale: 初始变形缩放系数。
    - max_scale: 变形缩放系数的最大值。
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')  # 创建 3D 坐标轴

    # 绘制原始结构
    def plot_original_structure():
        """绘制原始结构的节点和元素"""
        for i, element in enumerate(elements):
            # 获取节点坐标
            x_vals = [element.node1.x, element.node2.x]
            y_vals = [element.node1.y, element.node2.y]
            z_vals = [element.node1.z, element.node2.z]

            # 绘制元素
            ax.plot(x_vals, y_vals, z_vals, 'b-', linewidth=2, label="Original" if i == 0 else "")

            # 在元素中点处添加标签
            mid_x = (x_vals[0] + x_vals[1]) / 2
            mid_y = (y_vals[0] + y_vals[1]) / 2
            mid_z = (z_vals[0] + z_vals[1]) / 2
            ax.text(mid_x, mid_y, mid_z, f'({i + 1})', fontsize=10)

    plot_original_structure()

    # 设置变形后的结构的线条对象
    deformed_lines = [ax.plot([], [], [], 'r--', linewidth=2, label="Deformed" if i == 0 else "")[0]
                      for i in range(len(elements))]

    # 绘制节点
    def plot_nodes():
        """绘制节点"""
        for i, node in enumerate(nodes):
            ax.plot(node.x, node.y, node.z, 'ro')
            ax.text(node.x, node.y, node.z, f'{i + 1}', fontsize=12)

    plot_nodes()

    # 设置坐标轴范围并确保比例为1:1:1
    def set_axes_limits():
        """设置坐标轴范围并保持比例1:1:1"""
        node_x = [node.x for node in nodes]
        node_y = [node.y for node in nodes]
        node_z = [node.z for node in nodes]

        # 获取最小值和最大值
        x_min, x_max = min(node_x), max(node_x)
        y_min, y_max = min(node_y), max(node_y)
        z_min, z_max = min(node_z), max(node_z)

        # 计算坐标轴范围
        x_range = x_max - x_min
        y_range = y_max - y_min
        z_range = z_max - z_min
        max_range = max(x_range, y_range, z_range)

        # 设置坐标轴中心和边界
        x_mid = (x_max + x_min) / 2
        y_mid = (y_max + y_min) / 2
        z_mid = (z_max + z_min) / 2

        ax.set_xlim(x_mid - max_range / 2, x_mid + max_range / 2)
        ax.set_ylim(y_mid - max_range / 2, y_mid + max_range / 2)
        ax.set_zlim(z_mid - max_range / 2, z_mid + max_range / 2)

        ax.set_box_aspect([1, 1, 1])  # 保持1:1:1的比例

    set_axes_limits()

    # 设置图例和坐标轴标签
    ax.legend()
    ax.grid(True)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Original vs. Deformed System")

    # 创建缩放滑动条
    ax_scale = plt.axes((0.1, 0.1, 0.8, 0.03), facecolor='lightgoldenrodyellow')
    scale_slider = Slider(ax_scale, 'Scale', 0.1, max_scale, valinit=initial_scale, valstep=0.1)

    # 更新变形系统的绘制
    def update(val):
        """更新变形后的结构"""
        scale = scale_slider.val
        # 计算变形后的节点坐标
        deformed_x = [node.x + scale * U[n_dof_per_node * i] for i, node in enumerate(nodes)]
        deformed_y = [node.y + scale * U[n_dof_per_node * i + 1] for i, node in enumerate(nodes)]
        deformed_z = [node.z + scale * U[n_dof_per_node * i + 2] for i, node in enumerate(nodes)]

        # 更新每个元素的变形状态
        for i, element in enumerate(elements):
            x_values = [deformed_x[nodes.index(element.node1)], deformed_x[nodes.index(element.node2)]]
            y_values = [deformed_y[nodes.index(element.node1)], deformed_y[nodes.index(element.node2)]]
            z_values = [deformed_z[nodes.index(element.node1)], deformed_z[nodes.index(element.node2)]]
            deformed_lines[i].set_data(x_values, y_values)
            deformed_lines[i].set_3d_properties(z_values)

        fig.canvas.draw_idle()  # 更新图像

    # 将缩放滑动条与更新函数绑定
    scale_slider.on_changed(update)

    # 初次绘制变形后的结构
    update(initial_scale)

    plt.show()
