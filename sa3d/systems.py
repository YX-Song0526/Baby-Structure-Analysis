import numpy as np
from numpy import cos, sin
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from sa3d.elements import Node, Bar, BeamColumn
from sa3d.plot import plot_system


class Truss3D:
    def __init__(self):
        self.nodes: list[Node] = []
        self.bars: list[Bar] = []
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

    def cal_K_total(self):
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

    def solve_disp(self, tolerance=1e-10):
        n = len(self.nodes)
        free_dof = list((set(range(3 * n)).difference(self.fixed_dof)))

        K = self.cal_K_total()
        K_ff = K[np.ix_(free_dof, free_dof)]
        F_ff = np.array([self.F[i] for i in free_dof])
        U_ff = np.linalg.pinv(K_ff) @ F_ff

        U_ff[np.abs(U_ff) < tolerance] = 0

        U = np.zeros(3 * n)
        U[free_dof] = U_ff

        return U

    def plot(self, initial_scale=1.0, scale_max=10000.0):
        U = self.solve_disp()
        plot_system(self.bars,
                    self.nodes,
                    U,
                    3,
                    initial_scale,
                    scale_max)


class Frame3D:
    def __init__(self):
        self.nodes = []
        self.elements = []
        self.FnM = None
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
        self.FnM = []
        current_index = 0
        for node in self.nodes:
            self.node_indices.append(current_index)
            self.FnM += node.dof * [0.0]
            current_index += node.dof
        self.FnM = np.array(self.FnM)

    def add_single_force(self, node_idx: int, Fx=0., Fy=0., Fz=0.):
        if self.FnM is None:
            self.assign_dof()

        index = self.node_indices[node_idx - 1]
        self.FnM[index], self.FnM[index + 1], self.FnM[index + 2] = Fx, Fy, Fz

    def add_single_moment(self, node_idx: int, Mx=0., My=0., Mz=0.):

        if self.nodes[node_idx - 1].dof == 3:
            print("Warning: Can not add moment on a bar")
            return

        if self.FnM is None:
            self.assign_dof()

        index = self.node_indices[node_idx - 1]
        self.FnM[index + 3], self.FnM[index + 4], self.FnM[index + 5] = Mx, My, Mz

    def add_distributed_force(self, bc_idx: int, q: float):

        if bc_idx - 1 < 0 or bc_idx > len(self.elements) or (not isinstance(self.elements[bc_idx - 1], BeamColumn)):
            print(f"Error: BeamColumn({bc_idx}) does not exist.")
            return

        if self.FnM is None:
            self.assign_dof()

        bc = self.elements[bc_idx - 1]
        L, phi = bc.L, bc.Phi

        Fx_eq, Fy_eq, M_eq = -q * L * sin(phi) / 2.0, q * L * cos(phi) / 2.0, q * L ** 2 / 8.0

        index1, index2 = self.node_indices[self.nodes.index(bc.node1)], self.node_indices[self.nodes.index(bc.node2)]

        self.FnM[index1] += Fx_eq  # 起始节点的水平等效力
        self.FnM[index1 + 1] += Fy_eq  # 起始节点的竖直等效力
        self.FnM[index1 + 2] += M_eq  # 起始节点的等效弯矩

        self.FnM[index2] += Fx_eq  # 起始节点的水平等效力
        self.FnM[index2 + 1] += Fy_eq  # 起始节点的竖直等效力
        self.FnM[index2 + 2] -= M_eq  # 起始节点的等效弯矩

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

    def cal_K_total(self):

        n = len(self.FnM)
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

    def solve_disp(self):
        n = len(self.FnM)
        free_dof = list((set(range(n)).difference(self.fixed_dof)))

        K = self.cal_K_total()
        # print(K[np.ix_(free_dof, free_dof)])
        K_ff = csr_matrix(K[np.ix_(free_dof, free_dof)])
        F_ff = np.array([self.FnM[i] for i in free_dof])
        U_ff = spsolve(K_ff, F_ff)

        # U_ff[np.abs(U_ff) < tolerance] = 0

        U = np.zeros(n)
        U[free_dof] = U_ff

        return U

    def solve_reaction(self, tolerance=1e-8):
        n = len(self.FnM)
        K = csr_matrix(self.cal_K_total())
        U = self.solve_disp()
        Q = K @ U
        Q[np.abs(Q) < tolerance] = 0

        return Q

    def plot(self, initial_scale=1.0, scale_max=10000.0):
        U = self.solve_disp()
        plot_system(self.elements,
                    self.nodes,
                    U,
                    n_dof_per_node=6,
                    initial_scale=initial_scale,
                    max_scale=scale_max)
