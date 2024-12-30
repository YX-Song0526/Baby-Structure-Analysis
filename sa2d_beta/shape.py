from abc import ABC, abstractmethod
from numpy import pi


class Shape(ABC):
    def __init__(self):
        self.A = 0.0  # 面积
        self.I = 0.0  # 惯性矩
        self.y_max = 0.0  # y方向上的最大距离
        self.parameters = {}  # 存储形状参数

    @abstractmethod
    def update(self, **kwargs):
        """
        更新截面形状参数，采用关键字参数
        """
        pass  # 每个子类需要实现自己的更新方法

    @abstractmethod
    def recalculate(self):
        """
        根据当前参数重新计算面积、惯性矩和其他属性
        """
        pass

    @classmethod
    def get_parameters(cls) -> dict[str, str]:
        """
        返回每个形状类的参数名称和说明，方便用户了解哪些参数是可更新的。
        """
        raise NotImplementedError("子类需要实现 get_parameters 方法")


class Circle(Shape):
    def __init__(self, R: float):
        super().__init__()
        self.parameters['R'] = R
        self.update(R=R)

    def update(self, **kwargs):
        for key, value in kwargs.items():
            # 更新参数
            self.parameters[key] = value

        # 更新截面积，惯性矩和y_max
        self.recalculate()

    def recalculate(self):
        R = self.parameters['R']
        self.A = pi * R ** 2
        self.I = pi * R ** 4 / 4
        self.y_max = R

    @classmethod
    def get_parameters(cls) -> dict[str, str]:
        return {
            'R': '圆的半径'
        }


class Rectangle(Shape):
    def __init__(self, b: float, h: float):
        super().__init__()
        self.parameters['b'] = b
        self.parameters['h'] = h
        self.update(b=b, h=h)

    def update(self, **kwargs):
        for key, value in kwargs.items():
            # 更新参数
            self.parameters[key] = value
        # 更新面积、惯性矩和最大y值
        self.recalculate()

    def recalculate(self):
        b = self.parameters['b']
        h = self.parameters['h']
        self.A = b * h
        self.I = b * h ** 3 / 12
        self.y_max = h / 2

    @classmethod
    def get_parameters(cls) -> dict[str, str]:
        return {
            'b': '矩形的宽度',
            'h': '矩形的高度'
        }


class Box(Shape):
    def __init__(self,
                 a: float,
                 b: float,
                 t1: float,
                 t2: float,
                 t3: float,
                 t4: float):
        super().__init__()
        self.parameters['a'] = a
        self.parameters['b'] = b
        self.parameters['t1'] = t1
        self.parameters['t2'] = t2
        self.parameters['t3'] = t1
        self.parameters['t4'] = t2
        self.update(a=a, b=b, t1=t1, t2=t2, t3=t3, t4=t4)

    def update(self, **kwargs):
        for key, value in kwargs.items():
            # 更新参数
            self.parameters[key] = value
        # 更新面积、惯性矩和最大y值
        self.recalculate()

    def recalculate(self):
        a = self.parameters['a']
        b = self.parameters['b']
        t1 = self.parameters['t1']
        t2 = self.parameters['t2']
        t3 = self.parameters['t3']
        t4 = self.parameters['t4']
        self.A = a * b - (a - t1 - t3) * (b - t2 - t4)
        self.I = a * b ** 3 / 12 - (a - t1 - t3) * (b - t2 - t4) ** 3 / 12
        self.y_max = b / 2

    @classmethod
    def get_parameters(cls) -> dict[str, str]:
        return {
            'a': '宽度',
            'b': '高度',
            't1': '右壁厚',
            't2': '上壁厚',
            't3': '左壁厚',
            't4': '下壁厚'
        }


class Generalized(Shape):
    def __init__(self,
                 A: float,
                 I: float,
                 y_max: float):
        super().__init__()
        self.parameters['A'] = A
        self.parameters['I'] = I
        self.parameters['y_max'] = y_max
        self.update(A=A, I=I, y_max=y_max)

    def update(self, **kwargs):
        for key, value in kwargs.items():
            # 更新参数
            self.parameters[key] = value

        # 更新截面积，惯性矩和y_max
        self.recalculate()

    def recalculate(self):
        A = self.parameters['A']
        I = self.parameters['I']
        y_max = self.parameters['y_max']
        self.A = A
        self.I = I
        self.y_max = y_max
