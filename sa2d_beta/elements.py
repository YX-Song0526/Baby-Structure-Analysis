import numpy as np
from sa2d_beta.section import Section
from sa2d_beta.matrices import big_mat, trans_mat_for_frame


class Node:
    def __init__(self,
                 x: float,
                 y: float):
        self.x = x
        self.y = y


class Beam:
    def __init__(self,
                 node1: Node,
                 node2: Node,
                 section: Section):
        """
        二维钢架单元初始化

        Args:
            node1: 节点1
            node2: 节点2
            section: 截面

        """
        self.node1 = node1
        self.node2 = node2
        self.section = section

        # 计算单元长度
        self.L = np.sqrt((node2.x - node1.x) ** 2 +
                         (node2.y - node1.y) ** 2)

        # 偏转角（相对于整体x坐标的正方向）
        self.Phi = np.arctan2((node2.y - node1.y), (node2.x - node1.x))

        self.K_local = None
        self.K_global = None

        self.update()

    def update(self):
        E = self.section.material.E
        A = self.section.shape.A
        I = self.section.shape.I
        self.K_local = big_mat(E, A, I, self.L)
        self.K_global = trans_mat_for_frame(self.Phi).T @ self.K_local @ trans_mat_for_frame(self.Phi)

    def update_shape_params(self, **kwargs):
        self.section.update_shape_params(**kwargs)
        self.update()
