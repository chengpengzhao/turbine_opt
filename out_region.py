from __future__ import division
from math import *
import numpy as np
from ofblockmeshdicthelper import *
import matplotlib.pyplot as plt
from sympy import solve, symbols
import fileinput
from scipy.optimize import fsolve

# Python 3.6
# ========================================================================================
# mesh parameters
Rregion = 4
Length = 200
Width = 27
L0 = (Length - Width) / 2
#FCells = 70
#DCells = 100
#CellsCircle = 3 * (2 * FCells + 2 * DCells)
#CellsV = round(CellsCircle / 4)
CellsV=50
CellsP = 50
GradingP = 5
GradingV = 10
length_z = 0.1


# ========================================================================================
# calculate coordinates

def getInsertPointBetweenCircleAndLine(k, b, m, n, r):
    x = symbols('x')
    y = k * x + b
    Solve = solve((x - m) ** 2 + (y - n) ** 2 - r ** 2, x)
    return Solve


[x0, y0] = [-(L0 + Width / 2), -Width / 2]
[x7, y7] = [-(L0 + Width / 2), Width / 2]
[x1, y1] = [-(Width / 2), -Width / 2]
[x6, y6] = [-(Width / 2), Width / 2]
[x2, y2] = [(Width / 2), -Width / 2]
[x5, y5] = [(Width / 2), Width / 2]
[x3, y3] = [(L0 + Width / 2), -Width / 2]
[x4, y4] = [(L0 + Width / 2), Width / 2]
[x8, x10] = getInsertPointBetweenCircleAndLine(tan(pi / 4), -tan(pi / 4) * x5 + y5, 0, 0, Rregion)
[y8, y10] = [x8, x10]
[x11, x9] = getInsertPointBetweenCircleAndLine(-tan(pi / 4), tan(pi / 4) * x6 + y6, 0, 0, Rregion)
[y9, y11] = [-x9, -x11]
# test plot
# for i in range(0,12):
# plt.plot(locals()['x' + str(i)], locals()['y' + str(i)], 'ro', label="point")
# plt.gca().axis('equal')
# plt.show()
# ========================================================================================
# create blockmesh
# prepare ofblockmeshdicthelper.BlockMeshDict instance to
# gather vertices, blocks, faces and boundaries.
bmd = BlockMeshDict()
bmd.set_metric('cm')
basevs = [Vertex(x0, y0, 0, 'v0'),
          Vertex(x1, y1, 0, 'v1'),
          Vertex(x2, y2, 0, 'v2'),
          Vertex(x3, y3, 0, 'v3'),
          Vertex(x4, y4, 0, 'v4'),
          Vertex(x5, y5, 0, 'v5'),
          Vertex(x6, y6, 0, 'v6'),
          Vertex(x7, y7, 0, 'v7'),
          Vertex(x8, y8, 0, 'v8'),
          Vertex(x9, y9, 0, 'v9'),
          Vertex(x10, y10, 0, 'v10'),
          Vertex(x11, y11, 0, 'v11'),
          ]

for v in basevs:
    bmd.add_vertex(v.x, v.y, v.z, v.name + '-z')
    bmd.add_vertex(v.x, v.y, v.z + length_z, v.name + '+z')

# define circle and spline
z = length_z / 2
for i in ['+z', '-z']:
    # circle
    bmd.add_arcedge(('v11' + i, 'v10' + i), 'arc-1' + i, Vertex(0, Rregion, length_z / 2 + eval(i), 'arc-1' + i))
    bmd.add_arcedge(('v11' + i, 'v8' + i), 'arc-4' + i, Vertex(-Rregion, 0, length_z / 2 + eval(i), 'arc-4' + i))
    bmd.add_arcedge(('v9' + i, 'v10' + i), 'arc-2' + i, Vertex(Rregion, 0, length_z / 2 + eval(i), 'arc-2' + i))
    bmd.add_arcedge(('v8' + i, 'v9' + i), 'arc-3' + i, Vertex(0, -Rregion, length_z / 2 + eval(i), 'arc-3' + i))


def vnamegen(x0z0, x1z0, x1z1, x0z1):
    return (x0z0 + '-z', x1z0 + '-z', x1z1 + '-z', x0z1 + '-z',
            x0z0 + '+z', x1z0 + '+z', x1z1 + '+z', x0z1 + '+z')


# == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
# create blocks
b0 = bmd.add_hexblock(vnamegen('v0', 'v1', 'v6', 'v7'), (CellsP, CellsV, 1), 'b0', SimpleGrading(1 / GradingP, 1, 1))
b5 = bmd.add_hexblock(vnamegen('v2', 'v3', 'v4', 'v5'), (CellsP, CellsV, 1), 'b5', SimpleGrading(GradingP, 1, 1))
b4 = bmd.add_hexblock(vnamegen('v1', 'v2', 'v9', 'v8'), (CellsV, CellsV, 1), 'b4', SimpleGrading(1, 1 / GradingV, 1))
b3 = bmd.add_hexblock(vnamegen('v2', 'v5', 'v10', 'v9'), (CellsV, CellsV, 1), 'b3', SimpleGrading(1, 1 / GradingV, 1))
b2 = bmd.add_hexblock(vnamegen('v5', 'v6', 'v11', 'v10'), (CellsV, CellsV, 1), 'b2', SimpleGrading(1, 1 / GradingV, 1))
b1 = bmd.add_hexblock(vnamegen('v6', 'v1', 'v8', 'v11'), (CellsV, CellsV, 1), 'b1', SimpleGrading(1, 1 / GradingV, 1))
# == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
# boundary
bmd.add_boundary('patch', 'inlet', [b0.face('xm')])
bmd.add_boundary('patch', 'wall-l1', [b0.face('yp')])
bmd.add_boundary('patch', 'wall-h1', [b0.face('ym')])
bmd.add_boundary('patch', 'wall-l2', [b2.face('ym')])
bmd.add_boundary('patch', 'wall-h2', [b4.face('ym')])
bmd.add_boundary('patch', 'wall-l3', [b5.face('yp')])
bmd.add_boundary('patch', 'wall-h3', [b5.face('ym')])
bmd.add_boundary('patch', 'outlet', [b5.face('xp')])

bmd.add_boundary('patch', 'interface-1', [b4.face('yp')])
bmd.add_boundary('patch', 'interface-2', [b3.face('yp')])
bmd.add_boundary('patch', 'interface-3', [b2.face('yp')])
bmd.add_boundary('patch', 'interface-4', [b1.face('yp')])

# == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
# prepare for output
bmd.assign_vertexid()
f = open('outRegion_blockMeshDict', 'wb')
f.write(bmd.format().encode())
f.close()
