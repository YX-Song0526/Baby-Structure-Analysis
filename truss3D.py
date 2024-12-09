from sa3d.systems import Truss3D

s = Truss3D()

s.add_node(0, 0, 0)
s.add_node(2, 0, 0)
s.add_node(0, 2, 0)
s.add_node(2, 2, 0)
s.add_node(0, 0, 3)
s.add_node(2, 0, 3)
s.add_node(0, 2, 3)
s.add_node(2, 2, 3)

s.add_bar(1, 5)
s.add_bar(2, 6)
s.add_bar(3, 7)
s.add_bar(4, 8)
s.add_bar(5, 6)
s.add_bar(7, 8)
s.add_bar(5, 7)
s.add_bar(6, 8)
s.add_bar(1, 6)
s.add_bar(1,7)
s.add_bar(4, 6)
s.add_bar(4, 7)
# s.add_bar(5, 8)

s.add_fixed_sup(1, 2, 3, 4)
s.add_node_force(5, Fx=100000, Fy=100000)

print("节点位移：", s.solve_disp())
s.plot(scale_max=200)
