import numpy as np
from numpy import cos, sin
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from space.elements import Node, Bar, BeamColumn


class Truss3D:
    def __init__(self):
        self.nodes = []
        self.bars = []
        self.F = []
        self.fixed_dof = []

    def add_node(self,
                 x: float,
                 y: float,
                 z: float):

        self.nodes.append(Node(x, y, z))
        self.F += [0., 0., 0.]

    def add_bar(self,
                node1_index: int,
                node2_index: int,
                E=69e9,
                A=0.01):

        # 检查节点索引是否在范围内
        if node1_index - 1 < 0 or node2_index - 1 < 0 or node1_index > len(self.nodes) or node2_index > len(self.nodes):
            print(f"Error: One or both nodes for Bar({node1_index}, {node2_index}) do not exist.")
            return

        # 添加新的 Bar
        self.bars.append(Bar(self.nodes[node1_index - 1], self.nodes[node2_index - 1], E, A))

    def add_node_force(self,
                       node_index: int,
                       Fx=0.,
                       Fy=0.,
                       Fz=0.):

        self.F[3 * (node_index - 1)] = Fx
        self.F[3 * (node_index - 1) + 1] = Fy
        self.F[3 * (node_index - 1) + 2] = Fz

    def add_fixed_sup(self, *args):
        for i in args:
            self.fixed_dof += [3 * (i - 1), 3 * (i - 1) + 1, 3 * (i - 1) + 2]

    def cal_K(self):
        n = len(self.nodes)
        K = np.zeros((3 * n, 3 * n))
        for bar in self.bars:
            index1 = self.nodes.index(bar.node1)
            index2 = self.nodes.index(bar.node2)
            K_global = bar.K_global
            dof = [3 * index1, 3 * index1 + 1, 3 * index1 + 2, 3 * index2, 3 * index2 + 1, 3 * index2 + 2]
            for i in range(6):
                for j in range(6):
                    K[dof[i], dof[j]] += K_global[i, j]
        return K

    def solve(self, tolerance=1e-10):
        n = len(self.nodes)
        free_dof = list((set(range(3 * n)).difference(self.fixed_dof)))

        K = self.cal_K()
        K_ff = K[np.ix_(free_dof, free_dof)]
        F_ff = np.array([self.F[i] for i in free_dof])
        U_ff = np.linalg.pinv(K_ff) @ F_ff

        U_ff[np.abs(U_ff) < tolerance] = 0

        U = np.zeros(3 * n)
        U[free_dof] = U_ff

        return U

    def plot_system(self, initial_scale=1.0, scale_max=10000.0):
        # 计算初始位移向量
        U = self.solve()

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')  # 3D axes

        # 绘制原始系统
        for i, bar in enumerate(self.bars):
            x_values = [bar.node1.x, bar.node2.x]
            y_values = [bar.node1.y, bar.node2.y]
            z_values = [bar.node1.z, bar.node2.z]
            ax.plot(x_values, y_values, z_values, 'b-', linewidth=2, label="Original" if i == 0 else "")
            mid_x = (x_values[0] + x_values[1]) / 2
            mid_y = (y_values[0] + y_values[1]) / 2
            mid_z = (z_values[0] + z_values[1]) / 2
            ax.text(mid_x, mid_y, mid_z, f'({i + 1})', fontsize=10)

        # 设置变形后系统的绘制
        deformed_lines = []
        for i, bar in enumerate(self.bars):
            line, = ax.plot([], [], [], 'r--', linewidth=2, label="Deformed" if i == 0 else "")
            deformed_lines.append(line)

        # 绘制节点
        for i, node in enumerate(self.nodes):
            ax.plot(node.x, node.y, node.z, 'ro')
            ax.text(node.x, node.y, node.z, f'{i + 1}', fontsize=12)

        # 设置坐标范围（可以选择根据最大位移来动态计算）
        node_x = [node.x for node in self.nodes]
        node_y = [node.y for node in self.nodes]
        node_z = [node.z for node in self.nodes]
        x_min, x_max = min(node_x), max(node_x)
        y_min, y_max = min(node_y), max(node_y)
        z_min, z_max = min(node_z), max(node_z)

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_zlim(z_min, z_max)

        # 设置图例
        ax.legend()
        ax.grid(True)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title("Original vs. Deformed Truss System")

        # 创建滑动条，使用 scale_max 设置最大值
        ax_scale = plt.axes((0.1, 0.1, 0.8, 0.03), facecolor='lightgoldenrodyellow')
        scale_slider = Slider(ax_scale, 'Scale', 0.1, scale_max, valinit=initial_scale, valstep=0.1)

        # 更新变形结构的绘制
        def update(val):
            scale = scale_slider.val
            deformed_x = [node.x + scale * U[3 * i] for i, node in enumerate(self.nodes)]
            deformed_y = [node.y + scale * U[3 * i + 1] for i, node in enumerate(self.nodes)]
            deformed_z = [node.z + scale * U[3 * i + 2] for i, node in enumerate(self.nodes)]
            for i, bar in enumerate(self.bars):
                x_values = [deformed_x[self.nodes.index(bar.node1)], deformed_x[self.nodes.index(bar.node2)]]
                y_values = [deformed_y[self.nodes.index(bar.node1)], deformed_y[self.nodes.index(bar.node2)]]
                z_values = [deformed_z[self.nodes.index(bar.node1)], deformed_z[self.nodes.index(bar.node2)]]
                deformed_lines[i].set_data(x_values, y_values)
                deformed_lines[i].set_3d_properties(z_values)
            fig.canvas.draw_idle()  # 更新图像

        # 绑定滑动条到更新函数
        scale_slider.on_changed(update)

        # 初次绘制变形系统
        update(initial_scale)

        plt.show()


class Frame3D:
    def __init__(self):
        self.nodes = []
        self.elements = []
        self.QnM = None
        self.fixed_dof = []
        self.node_indices = []

    def add_node(self, x: float, y: float, z: float):
        self.nodes.append(Node(x, y, z))

    def add_bar(self, node1_idx: int, node2_idx: int, E, A):
        # 检查节点索引是否在范围内
        if node1_idx - 1 < 0 or node2_idx - 1 < 0 or node1_idx > len(self.nodes) or node2_idx > len(self.nodes):
            print(f"Error: One or both nodes for Bar({node1_idx}, {node2_idx}) do not exist.")
            return

        # 添加新的 bar
        self.elements.append(Bar(self.nodes[node1_idx - 1], self.nodes[node2_idx - 1], E, A))

    def add_beam_column(self, node1_idx: int, node2_idx: int, E: float, A: float, G: float, J: float, I: list):
        # 检查节点索引是否在范围内
        if node1_idx - 1 < 0 or node2_idx - 1 < 0 or node1_idx > len(self.nodes) or node2_idx > len(self.nodes):
            print(f"Error: One or both nodes for BeamColumn({node1_idx}, {node2_idx}) do not exist.")
            return

        self.nodes[node1_idx - 1].set_dof(6)
        self.nodes[node2_idx - 1].set_dof(6)

        # 添加新的beam column
        self.elements.append(BeamColumn(self.nodes[node1_idx - 1], self.nodes[node2_idx - 1], E, A, G, J, I))

    def assign_dof(self):
        self.QnM = []
        current_index = 0
        for node in self.nodes:
            self.node_indices.append(current_index)
            self.QnM += node.dof * [0.0]
            current_index += node.dof
        self.QnM = np.array(self.QnM)

    def add_single_force(self, node_idx: int, Fx=0., Fy=0., Fz=0.):
        if self.QnM is None:
            self.assign_dof()

        index = self.node_indices[node_idx - 1]
        self.QnM[index], self.QnM[index + 1], self.QnM[index + 2] = Fx, Fy, Fz

    def add_single_moment(self, node_idx: int, Mx=0., My=0., Mz=0.):

        if self.nodes[node_idx - 1].dof == 3:
            print("Warning: Can not add moment on a bar")
            return

        if self.QnM is None:
            self.assign_dof()

        index = self.node_indices[node_idx - 1]
        self.QnM[index + 3], self.QnM[index + 4], self.QnM[index + 5] = Mx, My, Mz

    def add_distributed_force(self, bc_idx: int, q: float):

        if bc_idx - 1 < 0 or bc_idx > len(self.elements) or (not isinstance(self.elements[bc_idx - 1], BeamColumn)):
            print(f"Error: BeamColumn({bc_idx}) does not exist.")
            return

        if self.QnM is None:
            self.assign_dof()

        bc = self.elements[bc_idx - 1]
        L, phi = bc.L, bc.Phi

        Fx_eq, Fy_eq, M_eq = -q * L * sin(phi) / 2.0, q * L * cos(phi) / 2.0, q * L ** 2 / 8.0

        index1, index2 = self.node_indices[self.nodes.index(bc.node1)], self.node_indices[self.nodes.index(bc.node2)]

        self.QnM[index1] += Fx_eq  # 起始节点的水平等效力
        self.QnM[index1 + 1] += Fy_eq  # 起始节点的竖直等效力
        self.QnM[index1 + 2] += M_eq  # 起始节点的等效弯矩

        self.QnM[index2] += Fx_eq  # 起始节点的水平等效力
        self.QnM[index2 + 1] += Fy_eq  # 起始节点的竖直等效力
        self.QnM[index2 + 2] -= M_eq  # 起始节点的等效弯矩

    def add_fixed_sup(self, *args):
        for i in args:
            index = sum(node.dof for node in self.nodes[:i - 1])
            if self.nodes[i - 1].dof == 3:
                self.fixed_dof += [index, index + 1, index + 2]
            elif self.nodes[i - 1].dof == 6:
                self.fixed_dof += [index, index + 1, index + 2,
                                   index + 3, index + 4, index + 5]

    def add_simple_sup(self, *args):
        for i in args:
            index = sum(node.dof for node in self.nodes[:i - 1])
            self.fixed_dof += [index, index + 1, index + 2]

    def cal_K(self):

        n = len(self.QnM)
        K = np.zeros((n, n))

        for element in self.elements:
            node1_idx, node2_idx = self.nodes.index(element.node1), self.nodes.index(element.node2)

            dof1_start, dof2_start = self.node_indices[node1_idx], self.node_indices[node2_idx]

            dof = []

            if isinstance(element, Bar):
                dof = list(range(dof1_start, dof1_start + 3)) + \
                      list(range(dof2_start, dof2_start + 3))

            elif isinstance(element, BeamColumn):
                dof = list(range(dof1_start, dof1_start + 6)) + \
                      list(range(dof2_start, dof2_start + 6))

            Ke = element.K_global  # 根据单元类型获取对应的刚度矩阵

            for i in range(len(dof)):
                for j in range(len(dof)):
                    K[dof[i], dof[j]] += Ke[i, j]

        return K

    def solve_una(self):
        n = len(self.QnM)
        free_dof = list((set(range(n)).difference(self.fixed_dof)))

        K = self.cal_K()
        # print(K[np.ix_(free_dof, free_dof)])
        K_ff = csr_matrix(K[np.ix_(free_dof, free_dof)])
        F_ff = np.array([self.QnM[i] for i in free_dof])
        U_ff = spsolve(K_ff, F_ff)

        # U_ff[np.abs(U_ff) < tolerance] = 0

        U = np.zeros(n)
        U[free_dof] = U_ff

        return U

    def solve_qnm(self, tolerance=1e-8):
        n = len(self.QnM)
        K = csr_matrix(self.cal_K())
        U = self.solve_una()
        Q = K @ U
        Q[np.abs(Q) < tolerance] = 0

        return Q

    def plot_system(self, initial_scale=1.0, scale_max=10000.0):
        # 计算初始位移向量
        U = self.solve_una()

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')  # 3D axes

        # 绘制原始系统
        for i, element in enumerate(self.elements):
            x_values = [element.node1.x, element.node2.x]
            y_values = [element.node1.y, element.node2.y]
            z_values = [element.node1.z, element.node2.z]
            ax.plot(x_values, y_values, z_values, 'b-', linewidth=2, label="Original" if i == 0 else "")
            mid_x = (x_values[0] + x_values[1]) / 2
            mid_y = (y_values[0] + y_values[1]) / 2
            mid_z = (z_values[0] + z_values[1]) / 2
            ax.text(mid_x, mid_y, mid_z, f'({i + 1})', fontsize=10)

        # 设置变形后系统的绘制
        deformed_lines = []
        for i, element in enumerate(self.elements):
            line, = ax.plot([], [], [], 'r--', linewidth=2, label="Deformed" if i == 0 else "")
            deformed_lines.append(line)

        # 绘制节点
        for i, node in enumerate(self.nodes):
            ax.plot(node.x, node.y, node.z, 'ro')
            ax.text(node.x, node.y, node.z, f'{i + 1}', fontsize=12)

        # 设置坐标范围，并确保比例为1:1:1，添加边距
        node_x = [node.x for node in self.nodes]
        node_y = [node.y for node in self.nodes]
        node_z = [node.z for node in self.nodes]
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
        ax.set_title("Original vs. Deformed Frame System")

        # 创建滑动条，使用 scale_max 设置最大值
        ax_scale = plt.axes((0.1, 0.1, 0.8, 0.03), facecolor='lightgoldenrodyellow')
        scale_slider = Slider(ax_scale, 'Scale', 0.1, scale_max, valinit=initial_scale, valstep=0.1)

        # 更新变形结构的绘制
        def update(val):
            scale = scale_slider.val
            deformed_x = [node.x + scale * U[6 * i] for i, node in enumerate(self.nodes)]
            deformed_y = [node.y + scale * U[6 * i + 1] for i, node in enumerate(self.nodes)]
            deformed_z = [node.z + scale * U[6 * i + 2] for i, node in enumerate(self.nodes)]
            for i, element in enumerate(self.elements):
                x_values = [deformed_x[self.nodes.index(element.node1)], deformed_x[self.nodes.index(element.node2)]]
                y_values = [deformed_y[self.nodes.index(element.node1)], deformed_y[self.nodes.index(element.node2)]]
                z_values = [deformed_z[self.nodes.index(element.node1)], deformed_z[self.nodes.index(element.node2)]]
                deformed_lines[i].set_data(x_values, y_values)
                deformed_lines[i].set_3d_properties(z_values)
            fig.canvas.draw_idle()  # 更新图像

        # 绑定滑动条到更新函数
        scale_slider.on_changed(update)

        # 初次绘制变形系统
        update(initial_scale)

        plt.show()
