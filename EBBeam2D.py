from structure2d import EB2D

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
s3.add_beam(1, 2, Iz=Iz)
s3.add_beam(2, 3, Iz=Iz)
s3.add_single_force(3, -3000)
s3.add_single_moment(2, 5000)
s3.add_single_moment(1, 5000)
s3.add_single_moment(3, 500)
s3.add_fixed_sup(1)
s3.add_simple_sup(2)

s3.plot_system(scale_max=200)
