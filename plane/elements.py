import numpy as np
from plane.matrices import little_mat, tiny_mat, big_mat, trans_mat_for_frame


class Node:
    def __init__(self, x: float, y: float, dof: int = 2):
        self.x = x
        self.y = y
        self.dof = dof

    def set_dof(self, new_dof: int):
        if new_dof < 0:
            print("Warning: The degree of freedom should be at least 0.")
            return
        self.dof = new_dof


class Bar:
    def __init__(self,
                 node1: Node,
                 node2: Node,
                 E: float,
                 A: float,
                 density: float = None):
        self.type = 'b'
        self.node1 = node1
        self.node2 = node2
        self.L = np.sqrt((node2.x - node1.x) ** 2 +
                         (node2.y - node1.y) ** 2)
        self.Phi = np.arctan2((node2.y - node1.y), (node2.x - node1.x))
        self.E = E
        self.A = A
        self.density = density
        self.K_local = self.E * self.A / self.L * np.array([[1, 0, -1, 0],
                                                            [0, 0, 0, 0],
                                                            [-1, 0, 1, 0],
                                                            [0, 0, 0, 0]])
        self.K_global = self.E * self.A / self.L * little_mat(self.Phi)
        self.Me = self.density * A * self.L / 6 * np.array([[2, 0, 1, 0],
                                                            [0, 2, 0, 1],
                                                            [1, 0, 2, 0],
                                                            [0, 1, 0, 2]]) if self.density is not None else None


class Beam:
    def __init__(self,
                 node1: Node,
                 node2: Node,
                 E: float,
                 Iz: float,
                 density: float = None):
        self.node1 = node1
        self.node2 = node2
        self.E = E
        self.Iz = Iz
        self.density = density
        self.L = node2.x - node1.x
        self.Ke = self.E * self.Iz / self.L ** 3 * tiny_mat(self.L)
        self.Me = self.density * self.A * self.L / 420 * \
                  np.array([[156, 22 * self.L, 54, -13 * self.L],
                            [22 * self.L, 4 * self.L ** 2, 13 * self.L, -3 * self.L ** 2],
                            [54, 13 * self.L, 156, -22 * self.L],
                            [-13 * self.L, -3 * self.L ** 2, -22 * self.L, 4 * self.L ** 2]]) \
            if self.density is not None else None


class BeamColumn:
    def __init__(self,
                 node1: Node,
                 node2: Node,
                 E: float,
                 A: float,
                 I: float,
                 density: float = None):
        self.type = 'bc'
        self.node1 = node1
        self.node2 = node2
        self.E = E
        self.A = A
        self.I = I
        self.density = density
        self.L = np.sqrt((node2.x - node1.x) ** 2 +
                         (node2.y - node1.y) ** 2)
        self.Phi = np.arctan2((node2.y - node1.y), (node2.x - node1.x))
        self.K_local = big_mat(self.E, self.A, self.I, self.L)
        self.K_global = trans_mat_for_frame(self.Phi).T @ self.K_local @ trans_mat_for_frame(self.Phi)
        self.Me = self.density * self.A * self.L * \
                  np.array([[1 / 3, 0, 0, 1 / 6, 0, 0],
                            [0, 13 / 35, 11 * self.L / 210, 0, 9 / 70, -13 * self.L / 420],
                            [0, 11 * self.L / 210, self.L ** 2 / 105, 0, 13 * self.L / 420, -self.L ** 2 / 140],
                            [1 / 6, 0, 0, 1 / 3, 0, 0],
                            [0, 9 / 70, 13 * self.L / 420, 0, 13 / 35, -11 * self.L / 210],
                            [0, -13 * self.L / 420, -self.L ** 2 / 140, 0, -11 * self.L / 210, self.L ** 2 / 105]]) \
            if self.density is not None else None
