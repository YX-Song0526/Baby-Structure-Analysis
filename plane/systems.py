import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.interpolate import CubicHermiteSpline
from typing import Union
from plane.components import Node, Bar, Beam, BeamColumn
from plane.matrix import trans_mat_for_bar, trans_mat_for_frame


class Truss2D:
    def __init__(self):
        self.nodes = []
        self.bars = []
        self.F = []
        self.fixed_dof = []

    def add_node(self,
                 x: float,
                 y: float):
        """添加节点"""
        self.nodes.append(Node(x, y))
        self.F += [0., 0.]

    def add_bar(self,
                node1_index: int,
                node2_index: int,
                E=69e9,
                A=0.01):
        """添加杆件"""
        
        # 检查节点索引是否在范围内
        if node1_index - 1 < 0 or node2_index - 1 < 0 or node1_index > len(self.nodes) or node2_index > len(self.nodes):
            print(f"Error: One or both nodes for Bar({node1_index}, {node2_index}) do not exist.")
            return

        # 添加新的 Bar
        self.bars.append(Bar(self.nodes[node1_index - 1], self.nodes[node2_index - 1], E, A))

    def add_node_force(self, node_index: int,
                       Fx=0.,
                       Fy=0.):
        """添加节点力"""
        self.F[2 * (node_index - 1)] = Fx
        self.F[2 * (node_index - 1) + 1] = Fy

    def add_fixed_sup(self, *args):
        """添加固定铰支座"""
        for i in args:
            self.fixed_dof += [2 * (i - 1), 2 * (i - 1) + 1]

    def add_movable_sup(self, *args, fixed: str):
        """添加活动铰支座"""
        for i in args:
            if fixed == 'x':
                self.fixed_dof += [2 * (i - 1)]
            elif fixed == 'y':
                self.fixed_dof += [2 * (i - 1) + 1]

    def cal_K_total(self):
        """计算系统总体刚度矩阵"""

        n = len(self.nodes)
        K = np.zeros((2 * n, 2 * n))
        for bar in self.bars:
            index1 = self.nodes.index(bar.node1)
            index2 = self.nodes.index(bar.node2)
            K_global = bar.K_global
            dof = [2 * index1, 2 * index1 + 1, 2 * index2, 2 * index2 + 1]
            for i in range(4):
                for j in range(4):
                    K[dof[i], dof[j]] += K_global[i, j]
        return K

    def solve(self, tolerance=1e-10):
        """求解节点位移"""

        n = len(self.nodes)
        free_dof = list((set(range(2 * n)).difference(self.fixed_dof)))

        K = self.cal_K_total()
        K_ff = K[np.ix_(free_dof, free_dof)]
        F_ff = np.array([self.F[i] for i in free_dof])
        U_ff = np.linalg.pinv(K_ff) @ F_ff

        U_ff[np.abs(U_ff) < tolerance] = 0

        U = np.zeros(2 * n)
        U[free_dof] = U_ff

        return U

    def plot_system(self,
                    initial_scale=1.0,
                    scale_max=10000.0):
        """可视化"""

        # 计算初始位移向量
        U = self.solve()

        fig, ax = plt.subplots()
        plt.subplots_adjust(left=0.1, bottom=0.25)  # 为滑动条留出空间

        # 绘制原始系统
        for i, bar in enumerate(self.bars):
            x_values = [bar.node1.x, bar.node2.x]
            y_values = [bar.node1.y, bar.node2.y]
            ax.plot(x_values, y_values, 'b-', linewidth=2, label="Original" if i == 0 else "")
            mid_x = (x_values[0] + x_values[1]) / 2
            mid_y = (y_values[0] + y_values[1]) / 2
            ax.text(mid_x, mid_y, f'({i + 1})', fontsize=10)

        # 设置变形后系统的绘制
        deformed_lines = []
        for i, bar in enumerate(self.bars):
            line, = ax.plot([], [], 'r--', linewidth=2, label="Deformed" if i == 0 else "")
            deformed_lines.append(line)

        # 绘制节点
        for i, node in enumerate(self.nodes):
            ax.plot(node.x, node.y, 'ro')
            ax.text(node.x, node.y, f'{i + 1}', fontsize=12)

        # 计算节点坐标范围
        node_x = [node.x for node in self.nodes]
        node_y = [node.y for node in self.nodes]
        x_min, x_max = min(node_x), max(node_x)
        y_min, y_max = min(node_y), max(node_y)

        # 添加一定的余量
        x_margin = (x_max - x_min) * 0.2  # 10% 的余量
        y_margin = (y_max - y_min) * 0.2
        ax.set_xlim(x_min - x_margin, x_max + x_margin)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)

        # 设置图例
        ax.legend()
        ax.grid(True)
        ax.set_aspect('equal')
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.title("Original vs. Deformed Truss System")

        # 创建滑动条，使用 scale_max 设置最大值
        ax_scale = plt.axes((0.1, 0.1, 0.8, 0.03), facecolor='lightgoldenrodyellow')
        scale_slider = Slider(ax_scale, 'Scale', 0.1, scale_max, valinit=initial_scale, valstep=0.1)

        # 更新变形结构的绘制
        def update(val):
            scale = scale_slider.val
            deformed_x = [node.x + scale * U[2 * i] for i, node in enumerate(self.nodes)]
            deformed_y = [node.y + scale * U[2 * i + 1] for i, node in enumerate(self.nodes)]
            for i, bar in enumerate(self.bars):
                x_values = [deformed_x[self.nodes.index(bar.node1)], deformed_x[self.nodes.index(bar.node2)]
                            ]
                y_values = [deformed_y[self.nodes.index(bar.node1)], deformed_y[self.nodes.index(bar.node2)]
                            ]
                deformed_lines[i].set_data(x_values, y_values)
            fig.canvas.draw_idle()  # 更新图像

        # 绑定滑动条到更新函数
        scale_slider.on_changed(update)

        # 初次绘制变形系统
        update(initial_scale)

        plt.show()


class EB2D:
    def __init__(self):
        self.nodes = []
        self.beams = []
        self.QnM = []
        self.fixed_dof = []

    def add_node(self, x: float, y: float):
        self.nodes.append(Node(x, y))
        self.QnM += [0., 0.]

    def add_beam(self, node1_index: int, node2_index: int, E=69e9, Iz=0.01):
        # 检查节点索引是否在范围内
        if node1_index - 1 < 0 or node2_index - 1 < 0 or node1_index > len(self.nodes) or node2_index > len(self.nodes):
            print(f"Error: One or both nodes for Bar({node1_index}, {node2_index}) do not exist.")
            return

        # 添加新的 Bar
        self.beams.append(Beam(self.nodes[node1_index - 1], self.nodes[node2_index - 1], E, Iz))

    def add_single_force(self, node_index: int, F=0.):
        self.QnM[2 * (node_index - 1)] = F

    def add_single_moment(self, node_index: int, M=0.):
        self.QnM[2 * node_index - 1] = M

    def add_distributed_force(self, beam_idx: int, q: float):

        if beam_idx - 1 < 0 or beam_idx > len(self.beams):
            print(f"Error: Beam({beam_idx}) does not exist.")
            return

        beam = self.beams[beam_idx - 1]
        L = beam.L

        F_eq = q * L / 2.0
        M_eq = q * L ** 2 / 8.0

        node1_idx = self.nodes.index(beam.node1)
        node2_idx = self.nodes.index(beam.node2)

        self.QnM[2 * node1_idx] += F_eq  # 起始节点的水平等效力
        self.QnM[2 * node2_idx] += F_eq  # 终止节点的水平等效力

        self.QnM[2 * node1_idx + 1] += M_eq  # 起始节点的等效弯矩
        self.QnM[2 * node2_idx + 1] -= M_eq  # 终止节点的等效弯矩

    def add_fixed_sup(self, *args):
        for i in args:
            self.fixed_dof += [2 * (i - 1), 2 * (i - 1) + 1]

    def add_simple_sup(self, *args):
        for i in args:
            self.fixed_dof += [2 * (i - 1)]

    def cal_K(self):
        n = len(self.nodes)
        K = np.zeros((2 * n, 2 * n))
        for beam in self.beams:
            index1 = self.nodes.index(beam.node1)
            index2 = self.nodes.index(beam.node2)
            Ke = beam.Ke
            dof = [2 * index1, 2 * index1 + 1, 2 * index2, 2 * index2 + 1]
            for i in range(4):
                for j in range(4):
                    K[dof[i], dof[j]] += Ke[i, j]
        return K

    def solve_una(self, tolerance=1e-10):
        n = len(self.nodes)
        free_dof = list((set(range(2 * n)).difference(self.fixed_dof)))

        K = self.cal_K()
        K_ff = K[np.ix_(free_dof, free_dof)]
        F_ff = np.array([self.QnM[i] for i in free_dof])
        U_ff = np.linalg.pinv(K_ff) @ F_ff

        U_ff[np.abs(U_ff) < tolerance] = 0

        U = np.zeros(2 * n)
        U[free_dof] = U_ff

        return U

    def solve_qnm(self, tolerance=1e-10):
        n = len(self.nodes)
        K = self.cal_K()
        U = self.solve_una()
        Q = K @ U
        Q[np.abs(Q) < tolerance] = 0

        return Q

    def plot_system(self, initial_scale=1.0, scale_max=1000.0):
        # 计算节点位移
        U = self.solve_una()

        # 创建图形和滑块
        fig, ax = plt.subplots()
        plt.subplots_adjust(left=0.1, bottom=0.3)  # 为滑动条留出空间

        # 绘制原始系统
        for i, beam in enumerate(self.beams):
            x_values = [beam.node1.x, beam.node2.x]
            y_values = [beam.node1.y, beam.node2.y]
            ax.plot(x_values, y_values, 'b-', linewidth=5, label="Original" if i == 0 else "")
            mid_x = (x_values[0] + x_values[1]) / 2
            mid_y = (y_values[0] + y_values[1]) / 2
            ax.text(mid_x, mid_y, f'({i + 1})', fontsize=10)

        # 绘制节点
        for i, node in enumerate(self.nodes):
            ax.plot(node.x, node.y, 'ro')
            ax.text(node.x, node.y, f'{i + 1}', fontsize=12)

        # 设置节点坐标范围
        node_x = [node.x for node in self.nodes]
        node_y = [node.y for node in self.nodes]
        x_min, x_max = min(node_x), max(node_x)
        y_min, y_max = min(node_y), max(node_y)

        # 设置坐标范围，并添加10%边距
        x_margin = (x_max - x_min) * 0.1 if x_max != x_min else 0.1
        y_margin = (y_max - y_min) * 0.5 if y_max != y_min else 0.5
        ax.set_xlim(x_min - x_margin, x_max + x_margin)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)

        # 设置图例、网格、比例和标题
        ax.grid(True)
        ax.set_aspect('equal')
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.title("Beam Deformation Visualization")

        # 创建滑动条，使用 scale_max 设置最大值
        ax_scale = plt.axes((0.1, 0.1, 0.8, 0.03), facecolor='lightgoldenrodyellow')
        scale_slider = Slider(ax_scale, 'Scale', 0.1, scale_max, valinit=initial_scale, valstep=0.1)

        # 初始化变形后系统的绘制，只创建一次
        deformed_lines = []
        for i, _ in enumerate(self.beams):
            if i == 0:
                line, = ax.plot([], [], 'r--', linewidth=2, label="Deformed")
            else:
                line, = ax.plot([], [], 'r--', linewidth=2)
            deformed_lines.append(line)

        # 设置图例一次，避免重复
        ax.legend(loc='best')

        def update(val):
            scale = scale_slider.val
            for i, beam in enumerate(self.beams):
                node1_idx = self.nodes.index(beam.node1)
                node2_idx = self.nodes.index(beam.node2)

                # 获取变形后节点的位置
                x_deformed = [beam.node1.x, beam.node2.x]
                y_deformed = [beam.node1.y + scale * U[2 * node1_idx],
                              beam.node2.y + scale * U[2 * node2_idx]]

                # 使用变形后的位置和转角来创建 Hermite 曲线
                dx_dy1 = U[2 * node1_idx + 1] * scale
                dx_dy2 = U[2 * node2_idx + 1] * scale

                # 定义 Hermite 插值对象
                hermite_spline = CubicHermiteSpline([x_deformed[0], x_deformed[1]],
                                                    [y_deformed[0], y_deformed[1]],
                                                    [dx_dy1, dx_dy2])

                # 生成插值点
                x_spline = np.linspace(x_deformed[0], x_deformed[1], 50)
                y_spline = hermite_spline(x_spline)

                # 更新线条数据而不更改标签
                deformed_lines[i].set_data(x_spline, y_spline)

            fig.canvas.draw_idle()  # 更新图像

        # 绑定滑动条到更新函数
        scale_slider.on_changed(update)

        # 初次绘制变形系统
        update(initial_scale)

        plt.show()


class Frame2D:
    def __init__(self):
        self.nodes: list[Node] = []
        self.elements: list[Union[Bar, BeamColumn]] = []
        self.Q = None
        self.fixed_dof = []
        self.node_indices = []

    def add_node(self,
                 x: float,
                 y: float):
        self.nodes.append(Node(x, y))

    def add_bar(self,
                node1_idx: int,
                node2_idx: int,
                E,
                A):
        # 检查节点索引是否在范围内
        if node1_idx - 1 < 0 or node2_idx - 1 < 0 or node1_idx > len(self.nodes) or node2_idx > len(self.nodes):
            print(f"Error: One or both nodes for Bar({node1_idx}, {node2_idx}) do not exist.")
            return

        # 添加新的 bar
        self.elements.append(Bar(self.nodes[node1_idx - 1], self.nodes[node2_idx - 1], E, A))

    def add_beam_column(self,
                        node1_idx: int,
                        node2_idx: int,
                        E,
                        A,
                        I):
        """节点1和节点2最好是从左向右布置"""

        # 检查节点索引是否在范围内
        if node1_idx - 1 < 0 or node2_idx - 1 < 0 or node1_idx > len(self.nodes) or node2_idx > len(self.nodes):
            print(f"Error: One or both nodes for BeamColumn({node1_idx}, {node2_idx}) do not exist.")
            return

        self.nodes[node1_idx - 1].set_dof(3)
        self.nodes[node2_idx - 1].set_dof(3)

        # 添加新的beam column
        self.elements.append(BeamColumn(self.nodes[node1_idx - 1], self.nodes[node2_idx - 1], E, A, I))

    def assign_dof(self):
        self.Q = []
        current_index = 0
        for node in self.nodes:
            self.node_indices.append(current_index)
            self.Q += node.dof * [0.0]
            current_index += node.dof
        self.Q = np.array(self.Q)

    def add_single_force(self, node_idx: int, Fx=0., Fy=0.):
        if self.Q is None:
            self.assign_dof()

        index = self.node_indices[node_idx - 1]
        self.Q[index], self.Q[index + 1] = Fx, Fy

    def add_single_moment(self, node_idx: int, M=0.):

        if self.nodes[node_idx - 1].dof == 2:
            print("Warning: Can not add moment on a bar")
            return

        if self.Q is None:
            self.assign_dof()

        index = self.node_indices[node_idx - 1]
        self.Q[index + 2] = M

    def add_distributed_force(self, bc_idx: int, q: float):

        if bc_idx - 1 < 0 or bc_idx > len(self.elements) or (not isinstance(self.elements[bc_idx - 1], BeamColumn)):
            print(f"Error: BeamColumn({bc_idx}) does not exist.")
            return

        if self.Q is None:
            self.assign_dof()

        bc = self.elements[bc_idx - 1]
        L, phi = bc.L, bc.Phi

        Q_eq_local = np.array([0, q * L / 2, q * L ** 2 / 8, 0, q * L / 2, -q * L ** 2 / 8])
        Q_eq_global = trans_mat_for_frame(phi) @ Q_eq_local

        print(Q_eq_global)

        index1, index2 = self.node_indices[self.nodes.index(bc.node1)], self.node_indices[self.nodes.index(bc.node2)]

        self.Q[index1] += Q_eq_global[0]  # 起始节点的水平等效力
        self.Q[index1 + 1] += Q_eq_global[1]  # 起始节点的竖直等效力
        self.Q[index1 + 2] += Q_eq_global[2]  # 起始节点的等效弯矩

        self.Q[index2] += Q_eq_global[3]  # 起始节点的水平等效力
        self.Q[index2 + 1] += Q_eq_global[4]  # 起始节点的竖直等效力
        self.Q[index2 + 2] += Q_eq_global[5]  # 起始节点的等效弯矩

    def add_fixed_sup(self, *args):
        for i in args:
            index = sum(node.dof for node in self.nodes[:i - 1])
            if self.nodes[i - 1].dof == 2:
                self.fixed_dof += [index, index + 1]
            elif self.nodes[i - 1].dof == 3:
                self.fixed_dof += [index, index + 1, index + 2]

    def add_simple_sup(self, *args):
        for i in args:
            index = sum(node.dof for node in self.nodes[:i - 1])
            self.fixed_dof += [index, index + 1]

    def cal_K(self):

        n = len(self.Q)
        K = np.zeros((n, n))

        for element in self.elements:
            node1_idx, node2_idx = self.nodes.index(element.node1), self.nodes.index(element.node2)

            dof1_start, dof2_start = self.node_indices[node1_idx], self.node_indices[node2_idx]

            dof = []

            if isinstance(element, Bar):
                dof = list(range(dof1_start, dof1_start + 2)) + \
                      list(range(dof2_start, dof2_start + 2))

            elif isinstance(element, BeamColumn):
                dof = list(range(dof1_start, dof1_start + 3)) + \
                      list(range(dof2_start, dof2_start + 3))

            Ke = element.K_global  # 根据单元类型获取对应的刚度矩阵

            for i in range(len(dof)):
                for j in range(len(dof)):
                    K[dof[i], dof[j]] += Ke[i, j]

        return K

    def solve_una(self, tolerance=1e-10):
        n = len(self.Q)
        free_dof = list((set(range(n)).difference(self.fixed_dof)))

        K = self.cal_K()

        K_ff = K[np.ix_(free_dof, free_dof)]
        F_ff = np.array([self.Q[i] for i in free_dof])
        U_ff = np.linalg.solve(K_ff, F_ff)

        U_ff[np.abs(U_ff) < tolerance] = 0

        U = np.zeros(n)
        U[free_dof] = U_ff

        return U

    def solve_qnm(self, tolerance=1e-10):
        n = len(self.Q)
        K = self.cal_K()
        U = self.solve_una()
        Q = K @ U
        Q[np.abs(Q) < tolerance] = 0

        return Q

    def get_dof_start(self, element: Union[Bar, BeamColumn]):
        # Extract the global indices for the degrees of freedom of this element's nodes
        node1_idx, node2_idx = self.nodes.index(element.node1), self.nodes.index(element.node2)
        dof1_start = self.node_indices[node1_idx]
        dof2_start = self.node_indices[node2_idx]

        return dof1_start, dof2_start

    def cal_element_nodal_force(self):
        U = self.solve_una()

        element_nodal_force_local = []

        for element in self.elements:

            dof1_start, dof2_start = self.get_dof_start(element)

            phi = element.Phi
            K_local = element.K_local

            if element.type == 'b':
                T_mat = trans_mat_for_bar(phi)
                u_global = U[dof1_start:dof1_start + 2].tolist() + U[dof2_start:dof2_start + 2].tolist()
                u_local = T_mat @ u_global
                f_local = K_local @ u_local

            elif element.type == 'bc':
                T_mat = trans_mat_for_frame(phi)
                u_global = U[dof1_start:dof1_start + 3].tolist() + U[dof2_start:dof2_start + 3].tolist()
                u_local = T_mat @ u_global
                f_local = K_local @ u_local

            else:
                f_local = None

            element_nodal_force_local.append(f_local)

        return element_nodal_force_local

    def get_element_deformed_pos(self, element: Union[Bar, BeamColumn], U: np.ndarray):
        n1_idx, n2_idx = self.nodes.index(element.node1), self.nodes.index(element.node2)
        index1, index2 = self.node_indices[n1_idx], self.node_indices[n2_idx]
        x = [element.node1.x + U[index1],
             element.node2.x + U[index2]]

        y = [element.node1.y + U[index1 + 1],
             element.node2.y + U[index2 + 1]]

        if element.type == 'b':
            return [x, y]
        elif element.type == 'bc':
            theta = [U[index1 + 2],
                     U[index2 + 2]]
            return [x, y, theta]

    def get_max_stress(self):
        R = 0.05
        ele_max_stress = []
        ele_nodal_force = self.cal_element_nodal_force()
        for i, f_e in enumerate(ele_nodal_force):
            dof = len(f_e)
            A = self.elements[i].A
            if dof == 4:
                Fx = f_e[0]
                axial_stress = abs(Fx / A)
                ele_max_stress.append(axial_stress)
            elif dof == 6:
                I = self.elements[i].I
                Fx = f_e[0]
                M1 = f_e[2]
                M2 = f_e[5]
                axial_stress = abs(Fx / A)
                bend_stress = max(abs(M1 * R / I), abs(M2 * R / I))
                stress = axial_stress + bend_stress
                ele_max_stress.append(stress)

        return ele_max_stress

    def plot_system(self, initial_scale=1.0, scale_max=1000.0):
        # 计算节点位移
        U = self.solve_una()

        fig, ax = plt.subplots(figsize=(8, 8))
        plt.subplots_adjust(left=0.1, bottom=0.3)  # 为滑动条留出空间

        # 绘制原始结构
        for i, element in enumerate(self.elements):
            x_values = [element.node1.x, element.node2.x]
            y_values = [element.node1.y, element.node2.y]
            color = 'b' if isinstance(element, Bar) else 'g'  # 区分Bar和BeamColumn
            ax.plot(x_values, y_values, color + '-', linewidth=5, label="Original" if i == 0 else "")
            mid_x = (x_values[0] + x_values[1]) / 2
            mid_y = (y_values[0] + y_values[1]) / 2
            ax.text(mid_x, mid_y, f'({i + 1})', fontsize=10)

        # 绘制节点
        for i, node in enumerate(self.nodes):
            ax.plot(node.x, node.y, 'ro')
            ax.text(node.x, node.y, f'{i + 1}', fontsize=12)

        # 设置节点坐标范围
        node_x = [node.x for node in self.nodes]
        node_y = [node.y for node in self.nodes]
        x_min, x_max = min(node_x), max(node_x)
        y_min, y_max = min(node_y), max(node_y)

        # 设置坐标范围，并添加10%边距
        x_margin = (x_max - x_min) * 0.3 if x_max != x_min else 0.3
        y_margin = (y_max - y_min) * 0.5 if y_max != y_min else 0.5
        ax.set_xlim(x_min - x_margin, x_max + x_margin)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)

        # 设置图例、网格、比例和标题
        ax.grid(True)
        ax.set_aspect('equal')
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.title("Frame System Deformation Visualization")

        # 创建滑动条，使用 scale_max 设置最大值
        ax_scale = plt.axes((0.1, 0.1, 0.8, 0.03), facecolor='lightgoldenrodyellow')
        scale_slider = Slider(ax_scale, 'Scale', 0.1, scale_max, valinit=initial_scale, valstep=0.1)

        # 初始化变形后系统的绘制，只创建一次
        deformed_lines = []
        for i, _ in enumerate(self.elements):
            if i == 0:
                line, = ax.plot([], [], 'r--', linewidth=2, label="Deformed")
            else:
                line, = ax.plot([], [], 'r--', linewidth=2)
            deformed_lines.append(line)

        # 设置图例一次，避免重复
        ax.legend(loc='best')

        def update(val):
            scale = scale_slider.val
            for i, element in enumerate(self.elements):

                node1_idx, node2_idx = self.nodes.index(element.node1), self.nodes.index(element.node2)

                index1, index2 = self.node_indices[node1_idx], self.node_indices[node2_idx]  # 获取自由度索引

                x_deformed = [element.node1.x + scale * U[index1],
                              element.node2.x + scale * U[index2]]

                y_deformed = [element.node1.y + scale * U[index1 + 1],
                              element.node2.y + scale * U[index2 + 1]]

                if isinstance(element, BeamColumn):

                    Phi = element.Phi

                    if Phi == 0:

                        # 对于水平的BeamColumn单元，使用变形后的位置和转角来创建 Hermite 曲线
                        dx_dy1 = U[index1 + 2] * scale
                        dx_dy2 = U[index2 + 2] * scale

                        hermite_spline = CubicHermiteSpline(
                            [x_deformed[0], x_deformed[1]],
                            [y_deformed[0], y_deformed[1]],
                            [dx_dy1, dx_dy2]
                        )

                        # 生成插值点
                        x_spline = np.linspace(x_deformed[0], x_deformed[1], 50)
                        y_spline = hermite_spline(x_spline)

                        deformed_lines[i].set_data(x_spline, y_spline)

                    else:

                        length = np.sqrt((x_deformed[1] - x_deformed[0]) ** 2 +
                                         (y_deformed[1] - y_deformed[0]) ** 2)
                        phi = np.arctan2(y_deformed[1] - y_deformed[0],
                                         x_deformed[1] - x_deformed[0])

                        dx_dy1 = U[index1 + 2] * scale + (Phi - phi)
                        dx_dy2 = U[index2 + 2] * scale + (Phi - phi)

                        # Hermite插值在局部坐标系中进行
                        hermite_spline = CubicHermiteSpline(
                            [0, length],
                            [0, 0],
                            [dx_dy1, dx_dy2]
                        )

                        # 局部坐标插值
                        x_old = np.linspace(0, length, 100)
                        y_old = hermite_spline(x_old)

                        # 坐标转换
                        x_new = x_deformed[0] + cos(phi) * x_old - sin(phi) * y_old
                        y_new = y_deformed[0] + sin(phi) * x_old + cos(phi) * y_old

                        deformed_lines[i].set_data(x_new, y_new)

                else:
                    # 其他单元绘制直线
                    deformed_lines[i].set_data(x_deformed, y_deformed)

            fig.canvas.draw_idle()

        # 绑定滑动条到更新函数
        scale_slider.on_changed(update)

        # 初次绘制变形系统
        update(initial_scale)

        plt.show()