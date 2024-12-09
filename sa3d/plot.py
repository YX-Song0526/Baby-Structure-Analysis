import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from sa3d.elements import Node, Bar, BeamColumn
from typing import Union


def plot_system(elements: list[Union[Bar, BeamColumn]],
                nodes: list[Node],
                U,
                n: int,
                initial_scale=1.0,
                scale_max=10000.0):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')  # 3D axes

    # 绘制原始系统
    for i, element in enumerate(elements):
        x_vals = [element.node1.x, element.node2.x]
        y_vals = [element.node1.y, element.node2.y]
        z_vals = [element.node1.z, element.node2.z]
        ax.plot(x_vals, y_vals, z_vals, 'b-', linewidth=2, label="Original" if i == 0 else "")
        mid_x = (x_vals[0] + x_vals[1]) / 2
        mid_y = (y_vals[0] + y_vals[1]) / 2
        mid_z = (z_vals[0] + z_vals[1]) / 2
        ax.text(mid_x, mid_y, mid_z, f'({i + 1})', fontsize=10)

    # 设置变形后系统的绘制
    deformed_lines = []
    for i, element in enumerate(elements):
        line, = ax.plot([], [], [], 'r--', linewidth=2, label="Deformed" if i == 0 else "")
        deformed_lines.append(line)

    # 绘制节点
    for i, node in enumerate(nodes):
        ax.plot(node.x, node.y, node.z, 'ro')
        ax.text(node.x, node.y, node.z, f'{i + 1}', fontsize=12)

    # 设置坐标范围，并确保比例为1:1:1，添加边距
    node_x = [node.x for node in nodes]
    node_y = [node.y for node in nodes]
    node_z = [node.z for node in nodes]
    x_min, x_max = min(node_x), max(node_x)
    y_min, y_max = min(node_y), max(node_y)
    z_min, z_max = min(node_z), max(node_z)

    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min
    max_range = max(x_range, y_range, z_range)

    x_mid = (x_max + x_min) / 2
    y_mid = (y_max + y_min) / 2
    z_mid = (z_max + z_min) / 2

    ax.set_xlim(x_mid - max_range / 2, x_mid + max_range / 2)
    ax.set_ylim(y_mid - max_range / 2, y_mid + max_range / 2)
    ax.set_zlim(z_mid - max_range / 2, z_mid + max_range / 2)

    ax.set_box_aspect([1, 1, 1])  # 确保坐标轴比例为 1:1:1

    # 设置图例
    ax.legend()
    ax.grid(True)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Original vs. Deformed System")

    # 创建滑动条，使用 scale_max 设置最大值
    ax_scale = plt.axes((0.1, 0.1, 0.8, 0.03), facecolor='lightgoldenrodyellow')
    scale_slider = Slider(ax_scale, 'Scale', 0.1, scale_max, valinit=initial_scale, valstep=0.1)

    # 更新变形结构的绘制
    def update(val):
        scale = scale_slider.val
        deformed_x = [node.x + scale * U[n * i] for i, node in enumerate(nodes)]
        deformed_y = [node.y + scale * U[n * i + 1] for i, node in enumerate(nodes)]
        deformed_z = [node.z + scale * U[n * i + 2] for i, node in enumerate(nodes)]
        for i, element in enumerate(elements):
            x_values = [deformed_x[nodes.index(element.node1)], deformed_x[nodes.index(element.node2)]]
            y_values = [deformed_y[nodes.index(element.node1)], deformed_y[nodes.index(element.node2)]]
            z_values = [deformed_z[nodes.index(element.node1)], deformed_z[nodes.index(element.node2)]]
            deformed_lines[i].set_data(x_values, y_values)
            deformed_lines[i].set_3d_properties(z_values)
        fig.canvas.draw_idle()  # 更新图像

    # 绑定滑动条到更新函数
    scale_slider.on_changed(update)

    # 初次绘制变形系统
    update(initial_scale)

    plt.show()
