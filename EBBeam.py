from sa2d.systems import EB2D

b = 0.1
h = 0.06
Iz = b * h ** 3 / 12

s1 = EB2D()
s1.add_node(0, 0)
s1.add_node(0.5, 0)
s1.add_node(1.0, 0)
s1.add_beam(1, 2, Iz=Iz)
s1.add_beam(2, 3, Iz=Iz)
# s1.add_single_force(2, F=-3000)
s1.add_distributed_force(2, -6000)
s1.add_distributed_force(1, -6000)
s1.add_simple_sup(1, 3)

s2 = EB2D()
s2.add_node(0., 0.)
s2.add_node(1.0, 0.)
s2.add_beam(1, 2, Iz=Iz)
s2.add_single_force(2, -2000)
s2.add_fixed_sup(1)

s3 = EB2D()
s3.add_node(0.0, 0.0)
s3.add_node(1.0, 0.0)
s3.add_node(2.0, 0.0)
s3.add_beam(1, 2)
s3.add_beam(2, 3)

s3.add_simple_sup(1, 3)
s3.add_distributed_force(1, -50000)
s3.add_distributed_force(2, -50000)


print("Displacement:", s3.solve_disp())
print("Reaction:", s3.solve_reaction())
s3.plot_system(scale_max=10000)
