from sa3d.systems import Truss3D

s = Truss3D()

s.add_node(0, 0, 0)
s.add_node(1, 0, 0)
s.add_node(0, 1, 0)
s.add_node(0, 0, 1)


s.add_bar(1, 2, E=69e9, A=0.02)
s.add_bar(1, 3, E=69e9, A=0.02)
s.add_bar(1, 4, E=69e9, A=0.02)
s.add_bar(2, 3, E=69e9, A=0.02)
s.add_bar(2, 4, E=69e9, A=0.02)
s.add_bar(3, 4, E=69e9, A=0.02)

s.add_fixed_sup(1, 2, 3)

s.add_node_force(4, Fz=50000)

print("节点位移：", s.solve())
s.plot_system(scale_max=5000)
