from plane.elements import *

n1 = Node(0, 0)
n2 = Node(4, 3)

e1 = BeamColumn(n1, n2, 1, 1, 1, 1)
e2 = Bar(n1, n2, 1, 1, 1)

print(e1.Me)