import numpy as np
from sa3d.matrices import augment, T_mat


class Node:
    def __init__(self,
                 x: float,
                 y: float,
                 z: float = 0,
                 dof: int = 3):
        self.x = x
        self.y = y
        self.z = z
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
                 E: float = 1.,
                 A: float = 1.):
        self.node1 = node1
        self.node2 = node2
        self.L = np.sqrt((node2.x - node1.x) ** 2 +
                         (node2.y - node1.y) ** 2 +
                         (node2.z - node1.z) ** 2)

        self.T_mat = np.array(
            [[(node2.x - node1.x) / self.L, (node2.y - node1.y) / self.L, (node2.z - node1.z) / self.L, 0, 0, 0],
             [0, 0, 0, (node2.x - node1.x) / self.L, (node2.y - node1.y) / self.L, (node2.z - node1.z) / self.L]]
        )
        self.E = E
        self.A = A
        self.K_local = self.E * self.A / self.L * np.array([[1, -1],
                                                            [-1, 1]])
        self.K_global = self.T_mat.T @ self.K_local @ self.T_mat


class BeamColumn:
    def __init__(self,
                 node1: Node,
                 node2: Node,
                 E: float,
                 A: float,
                 G: float,
                 J: float,
                 I: list):

        self.node1 = node1
        self.node2 = node2
        self.E = E
        self.A = A
        self.I = I
        self.G = G
        self.J = J

        dx = node2.x - node1.x
        dy = node2.y - node1.y
        dz = node2.z - node1.z

        T_wave = T_mat(dx, dy, dz)

        self.L = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        self.K_local = augment(self.E, self.A, self.G, self.J, self.I, self.L)
        self.K_global = T_wave.T @ self.K_local @ T_wave
