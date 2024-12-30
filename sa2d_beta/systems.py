import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.interpolate import CubicHermiteSpline
from sa2d_beta.elements import Node, Beam
from sa2d_beta.matrices import trans_mat_for_frame
from sa2d_beta.section import Section


class Frame2D:
    def __init__(self):
        self.nodes: list[Node] = []
        self.elements: list[Beam] = []
        self.FnM: list[float] = []
        self.fixed_dof: list[int] = []
        # self.node_indices = []

    def add_node(self,
                 x: float,
                 y: float):
        """
        添加节点

        Args:
            x: x 坐标
            y: y 坐标

        """
        self.nodes.append(Node(x, y))
        self.FnM += [0.0, 0.0, 0.0]  # 给节点力向量分配自由度

    def add_element(self,
                    node1_id: int,
                    node2_id: int,
                    section: Section):
        """
        添加单元

        Args:
            node1_id: 节点 1 编号
            node2_id: 节点 2 编号
            section: 单元截面

        """

        # 检查节点索引是否在范围内
        if node1_id - 1 < 0 or node2_id - 1 < 0 or node1_id > len(self.nodes) or node2_id > len(self.nodes):
            print(f"Error: One or both nodes for BeamColumn({node1_id}, {node2_id}) do not exist.")
            return

        beam = Beam(self.nodes[node1_id - 1],
                    self.nodes[node2_id - 1],
                    section)

        # 添加新的梁
        self.elements.append(beam)

    def add_single_force(self,
                         node_id: int,
                         Fx=0.0,
                         Fy=0.0):
        """
        添加节点力
        Args:
            node_id: 节点编号
            Fx: x 方向力
            Fy: y 方向力

        """

        i = node_id - 1
        self.FnM[3 * i], self.FnM[3 * i + 1] = Fx, Fy

    def add_single_moment(self, node_id: int, M=0.0):
        """
        添加集中弯矩

        Args:
            node_id: 节点编号
            M: 弯矩，正负表示方向

        """

        i = node_id - 1
        self.FnM[3 * i + 2] = M

    def add_fixed_sup(self, *args):
        """添加固定支座"""
        for node_id in args:
            i = node_id - 1
            self.fixed_dof += [3 * i, 3 * i + 1, 3 * i + 2]

    def add_simple_sup(self, *args):
        """添加简单支座"""
        for node_id in args:
            i = node_id - 1
            self.fixed_dof += [3 * i, 3 * i + 1]

    def cal_K_total(self):
        """计算总体刚度矩阵"""

        # 初始化总刚
        n = len(self.FnM)
        K = np.zeros((n, n))

        for element in self.elements:
            # 获取两节点的起始自由度编号
            i, j = self.nodes.index(element.node1), self.nodes.index(element.node2)

            # 自由度索引
            dof = [3 * i, 3 * i + 1, 3 * i + 2, 3 * j, 3 * j + 1, 3 * j + 2]

            Ke = element.K_global  # 单刚

            for i in range(6):
                for j in range(6):
                    K[dof[i], dof[j]] += Ke[i, j]

        return K

    def solve_disp(self, tolerance=1e-10):
        """
        求解节点位移

        Args:
            tolerance: 小于该值的位移将会被认为是0

        Returns:

        """
        n = len(self.FnM)
        free_dof = list((set(range(n)).difference(self.fixed_dof)))

        K = self.cal_K_total()

        K_ff = K[np.ix_(free_dof, free_dof)]
        F_f = np.array([self.FnM[i] for i in free_dof])
        U_f = np.linalg.solve(K_ff, F_f)

        U_f[np.abs(U_f) < tolerance] = 0

        U = np.zeros(n)
        U[free_dof] = U_f

        return U

    def solve_reaction(self, tolerance=1e-10):
        """
        求解反力

        Args:
            tolerance: 小于该值的力将会被认为是0

        Returns:

        """
        K = self.cal_K_total()
        U = self.solve_disp()
        Q = K @ U
        f = np.array(self.FnM)
        R = Q - f
        R[np.abs(R) < tolerance] = 0

        return R

    def get_element_dof(self, element):
        i, j = self.nodes.index(element.node1), self.nodes.index(element.node2)
        element_dof = [3 * i, 3 * i + 1, 3 * i + 2, 3 * j, 3 * j + 1, 3 * j + 2]
        return element_dof

    def cal_element_nodal_force(self):
        """
        求解单元的节点力
        """
        U = self.solve_disp()

        element_nodal_force_local = []

        for element in self.elements:
            # 获取单元的自由度索引
            dof_ids = self.get_element_dof(element)

            phi = element.Phi
            K_local = element.K_local

            # 转换矩阵
            T_mat = trans_mat_for_frame(phi)

            # 全局坐标系下的单元节点位移解
            u_global = U[dof_ids]

            # 局部坐标下的单元节点位移解
            u_local = T_mat @ u_global

            # 局部坐标系下的单元节点力
            f_local = K_local @ u_local

            element_nodal_force_local.append(f_local)

        return element_nodal_force_local

    def get_max_stress(self):
        """
        求解单元最大应力数组
        """
        ele_max_stress = []
        ele_nodal_force = self.cal_element_nodal_force()
        for i, f_e in enumerate(ele_nodal_force):
            # 获取单元的截面信息
            A = self.elements[i].section.shape.A
            I = self.elements[i].section.shape.I
            y_max = self.elements[i].section.shape.y_max

            # 获取计算应力需要用到的数据Fx, M1, M2
            Fx = f_e[0]
            M1 = f_e[2]
            M2 = f_e[5]

            # 拉压应力
            axial_stress = abs(Fx / A)

            # 弯曲应力
            bend_stress = max(abs(M1 * y_max / I), abs(M2 * y_max / I))

            # 最大总应力
            stress = axial_stress + bend_stress
            ele_max_stress.append(stress)

        return max(ele_max_stress)

    def plot_system(self, initial_scale=1.0, scale_max=1000.0):
        # 计算节点位移
        U = self.solve_disp()

        fig, ax = plt.subplots(figsize=(8, 8))
        plt.subplots_adjust(left=0.1, bottom=0.3)  # 为滑动条留出空间

        # 绘制原始结构
        for i, element in enumerate(self.elements):
            x_values = [element.node1.x, element.node2.x]
            y_values = [element.node1.y, element.node2.y]
            color = 'lightblue'
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

                if isinstance(element, Beam):

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
