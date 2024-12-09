from sa2d.systems import Frame2D

# -----------------------------------------------------------
#                         System1
# -----------------------------------------------------------

s1 = Frame2D()

s1.add_node(0., 0.)
s1.add_node(3., 0.)
s1.add_node(3., 3.)

s1.add_bar(1, 3, E=210e9, A=0.001)

s1.add_beam_column(1, 2, E=210e9, A=0.002, I=5e-5)

s1.add_fixed_sup(2, 3)
s1.add_single_force(1, Fy=-500000)

# -----------------------------------------------------------
#                         System2
# -----------------------------------------------------------

# s2 = Frame2D()
# s2.add_node(-1.0, 1.0)
# s2.add_node(0.0, 1.0)
# s2.add_node(1.0, 1.0)
# s2.add_node(0.0, 0.0)
# s2.add_beam_column(1, 2, E=210e9, A=20e-4, I=300e-8)
# s2.add_beam_column(3, 2, E=210e9, A=20e-4, I=300e-8)
# s2.add_beam_column(2, 4, E=210e9, A=20e-4, I=300e-8)
# s2.add_fixed_sup(1, 4)
# s2.add_simple_sup(3)
# s2.add_single_force(2, Fx=20000, Fy=-400000)
# s2.add_single_moment(2, -40000)

# -----------------------------------------------------------
#                         System3
# -----------------------------------------------------------

s3 = Frame2D()
s3.add_node(0.0, 8.0)
s3.add_node(12.0, 8.0)
s3.add_node(0.0, 0.0)
s3.add_node(12.0, 0.0)
s3.add_beam_column(1, 2, E=30e6, A=6.8, I=65)
s3.add_beam_column(3, 1, E=30e6, A=6.8, I=65)
s3.add_beam_column(2, 4, E=30e6, A=6.8, I=65)
s3.add_single_force(1, Fx=3000)

s3.add_distributed_force(1, -500.0)

s3.add_fixed_sup(3, 4)

# -----------------------------------------------------------
#                         System4
# -----------------------------------------------------------


# s4 = Frame2D()
# node_list = [(0., 8.), (4., 8.), (8., 8.), (2., 6.),
#              (4., 6.), (6., 6.), (0., 4.), (3., 4.),
#              (5., 4.), (8., 4.), (2., 2.), (4., 2.),
#              (6., 2.), (0., 0.), (4., 0.), (8., 0.)]
# element_list = [(1, 2), (2, 3), (7, 1), (5, 2), (10, 3),
#                 (8, 4), (9, 6), (7, 8),
#                 (9, 10), (14, 7), (11, 8), (13, 9), (16, 10),
#                 (15, 12), (14, 15), (15, 16)]
#
# for x, y in node_list:
#     s4.add_node(x, y)
#
# s4.add_node(0., 9.)
# s4.add_node(4., 9.)
# s4.add_node(8., 9.)
#
# for i, j in element_list:
#     s4.add_beam_column(i, j, 69e9, 0.01, 5e-5)
#
# s4.add_beam_column(4, 5, 2100e9, 0.01, 5e-5)
# s4.add_beam_column(5, 6, 2100e9, 0.01, 5e-5)
# s4.add_beam_column(11, 12, 2100e9, 0.01, 5e-5)
# s4.add_beam_column(12, 13, 2100e9, 0.01, 5e-5)
# s4.add_bar(17, 18, 69e9, 0.01)
# s4.add_bar(18, 19, 69e9, 0.01)
# s4.add_bar(1, 17, 69e9, 0.01)
# s4.add_bar(2, 18, 69e9, 0.01)
# s4.add_bar(3, 19, 69e9, 0.01)
# s4.add_bar(1, 18, 69e9, 0.01)
# s4.add_bar(17, 2, 69e9, 0.01)
# s4.add_bar(2, 19, 69e9, 0.01)
# s4.add_bar(18, 3, 69e9, 0.01)
#
# s4.add_single_force(17, Fy=-50000)
# s4.add_single_force(18, Fy=-8000)
# s4.add_single_force(19, Fy=-50000)
# s4.add_simple_sup(14, 15, 16)

# -----------------------------------------------------------
#                     Solution OUTPUT
# -----------------------------------------------------------

print("节点位移：", s1.solve_disp())
# print("节点力：", s3.solve_qnm())
#
# local_force = s3.cal_element_nodal_force()
# max_stress = s3.get_max_stress()
# print('单元节点力（局部）:', local_force)
# print('各单元最大应力:', max_stress)
s1.plot_system(scale_max=100)

# print(s3.elements[0].K_local)

