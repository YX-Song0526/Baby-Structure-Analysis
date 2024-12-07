from plane.systems import Truss2D

s1 = Truss2D()
s2 = Truss2D()

node_list = [(-0.5, 0.), (0.5, 0), (0., -1.), (0., -1.5)]
bar_list = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6),
            (1, 7), (1, 8), (2, 8), (2, 9), (3, 9), (3, 10), (4, 9), (4, 10), (5, 10), (5, 11), (6, 11), (6, 12),
            (7, 8), (8, 9), (9, 10), (10, 11), (11, 12)]

for i in range(6):
    s1.add_node(3 * i, 4)
for i in range(6):
    s1.add_node(3 * i, 0.)

for i, j in bar_list:
    s1.add_bar(i, j)

for x, y in node_list:
    s2.add_node(x, y)

s2.add_bar(1, 3)
s2.add_bar(3, 2)
s2.add_bar(3, 4)

s1.add_node_force(9, Fy=-7000)
s1.add_node_force(11, Fy=-6000)
s1.add_node_force(6, Fx=-8000)
s1.add_fixed_sup(7)
s1.add_movable_sup(12, fixed='y')

s2.add_node_force(4, Fy=-10000, Fx=500000)
s2.add_node_force(3, Fx=5000)
s2.add_fixed_sup(1)
s2.add_movable_sup(2, fixed='y')
U = s1.solve_disp()
U2 = s2.solve_disp()

# print("刚度矩阵:\n", s2.cal_K())
print("节点位移:\n", U)

s1.plot_system(scale_max=5000)
