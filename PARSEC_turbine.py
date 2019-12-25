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
def pcoef(
        xte, yte, rle,
        x_cre, y_cre, d2ydx2_cre, th_cre,
        surface):
    """evaluate the PARSEC coefficients"""
    # Initialize coefficients
    coef = np.zeros(6)

    # 1st coefficient depends on surface (pressure or suction)
    if surface.startswith('p'):
        coef[0] = -sqrt(2 * rle)
    else:
        coef[0] = sqrt(2 * rle)

    # Form system of equations
    A = np.array([
        [xte ** 1.5, xte ** 2.5, xte ** 3.5, xte ** 4.5, xte ** 5.5],
        [x_cre ** 1.5, x_cre ** 2.5, x_cre ** 3.5, x_cre ** 4.5,
         x_cre ** 5.5],
        [1.5 * sqrt(xte), 2.5 * xte ** 1.5, 3.5 * xte ** 2.5,
         4.5 * xte ** 3.5, 5.5 * xte ** 4.5],
        [1.5 * sqrt(x_cre), 2.5 * x_cre ** 1.5, 3.5 * x_cre ** 2.5,
         4.5 * x_cre ** 3.5, 5.5 * x_cre ** 4.5],
        [0.75 * (1 / sqrt(x_cre)), 3.75 * sqrt(x_cre), 8.75 * x_cre ** 1.5,
         15.75 * x_cre ** 2.5, 24.75 * x_cre ** 3.5]
    ])

    B = np.array([
        [yte - coef[0] * sqrt(xte)],
        [y_cre - coef[0] * sqrt(x_cre)],
        [tan(th_cre * pi / 180) - 0.5 * coef[0] * (1 / sqrt(xte))],
        [-0.5 * coef[0] * (1 / sqrt(x_cre))],
        [d2ydx2_cre + 0.25 * coef[0] * x_cre ** (-1.5)]
    ])

    # Solve system of linear equations
    X = np.linalg.solve(A, B)

    # Gather all coefficients
    coef[1:6] = X[0:5, 0]

    # Return coefficients
    return coef


# ========================================================================================
# Sample coefficients  NACA0012
'''
leading edge radius (r_le)
pressure and suction surface crest locations (x_pre, y_pre, x_suc, y_suc)
curvatures at the pressure and suction surface crest locations (d2y/dx2_pre, d2y/dx2_suc)
trailing edge coordinates (y_te)
trailing edge angles between the pressure and suction surface and the horizontal axis (th_pre, th_suc)
'''
rle_suc = .014927
rle_pre = .014181
x_suc = .29866
y_suc = .059404
x_pre = .29962
y_pre = -0.059632
d2ydx2_suc = -.42399
d2ydx2_pre = 0.445281
yte = 0
th_suc = -7.672047
th_pre = 7.59506

xte = 1.0  # Trailing edge x position
# Evaluate pressure (lower) surface coefficients
cf_pre = pcoef(xte, yte, rle_pre,
               x_pre, y_pre, d2ydx2_pre, th_pre,
               'pre')
# Evaluate suction (upper) surface coefficients
cf_suc = pcoef(xte, yte, rle_suc,
               x_suc, y_suc, d2ydx2_suc, th_suc,
               'suc')
x = (1 - np.cos(np.linspace(0, 1, int(np.ceil(1e3))) * np.pi)) / 2
uppery = np.array([0] * len(x))
lowery = np.array([0] * len(x))
for i in range(1, 7):
    uppery = uppery + cf_suc[i - 1] * x ** (i - 1 / 2)
    lowery = lowery + cf_pre[i - 1] * x ** (i - 1 / 2)

# ========================================================================================
ChordL = 2.54  # cm
x = ChordL * x;
uppery = ChordL * uppery;
lowery = ChordL * lowery;
# Plot this airfoil
# plt.plot(x, uppery, 'o--', markerfacecolor='red')
# plt.plot(x, lowery, 'o--', markerfacecolor='blue')
# plt.gca().axis('equal')
# plt.show()
# ========================================================================================
# mesh parameters
Rregion = 4
Rairfoil = 3.414
Rinner = 2.4
Lsquare = 2
beta = 15 * pi / 180  # pitch angle
theta = 120 * pi / 180
length_z = 0.1
RendTostart = 1.2 * Rregion
leadregionNumber = round(0.25 / 1 * len(x))

# mesh parameters
FCells = 70
DCells = 100
CCells = 70
UpCells = 70
DGrading = 10
UpGrading = 10
LeadGrading = 10
CGrading = 3
# ofblockmeshdicthelper doesn't support for this format
# ( length% cells% expansion )
FGrading = '((0.3 0.4 5)(0.4 0.2 1)(0.3 0.4 0.2))'  # (10*0.1=1) must be 1

# ========================================================================================
# calculate coordinates

# translation
[blade1upperx, blade1uppery] = [x - 0.25 * ChordL, uppery]
[blade1lowerx, blade1lowery] = [x - 0.25 * ChordL, lowery]

# rotation&translation
[blade1upperx, blade1uppery] = [blade1upperx * cos(-beta) - blade1uppery * sin(-beta),
                                blade1upperx * sin(-beta) + blade1uppery * cos(-beta) + Rairfoil]
[blade1lowerx, blade1lowery] = [blade1lowerx * cos(-beta) - blade1lowery * sin(-beta),
                                blade1lowerx * sin(-beta) + blade1lowery * cos(-beta) + Rairfoil]


def getInsertPointBetweenCircleAndLine(k, b, m, n, r):
    x = symbols('x')
    y = k * x + b
    Solve = solve((x - m) ** 2 + (y - n) ** 2 - r ** 2, x)
    return Solve


# ========================================================================================
# define points
[x5, y5] = [blade1upperx[0], blade1uppery[0]]
[x2, y2] = [blade1upperx[-1], blade1uppery[-1]]
x0 = max(getInsertPointBetweenCircleAndLine(tan(pi / 2 - beta), -tan(pi / 2 - beta) * x5 + y5, 0, 0, Rinner))
y0 = tan(pi / 2 - beta) * (x0 - x5) + y5
x4 = max(getInsertPointBetweenCircleAndLine(tan(pi / 2 - beta), -tan(pi / 2 - beta) * x5 + y5, 0, 0, Rregion))
y4 = tan(pi / 2 - beta) * (x4 - x5) + y5

x1 = max(getInsertPointBetweenCircleAndLine(tan(pi / 2 - beta), -tan(pi / 2 - beta) * x2 + y2, 0, 0, Rinner))
y1 = tan(pi / 2 - beta) * (x1 - x2) + y2
x3 = max(getInsertPointBetweenCircleAndLine(tan(pi / 2 - beta), -tan(pi / 2 - beta) * x2 + y2, 0, 0, Rregion))
y3 = tan(pi / 2 - beta) * (x3 - x2) + y2

[blade2lowerx, blade2lowery] = [blade1lowerx * cos(-theta) - blade1lowery * sin(-theta),
                                blade1lowerx * sin(-theta) + blade1lowery * cos(-theta)]
[blade2upperx, blade2uppery] = [blade1upperx * cos(-theta) - blade1uppery * sin(-theta),
                                blade1upperx * sin(-theta) + blade1uppery * cos(-theta)]
[blade3lowerx, blade3lowery] = [blade1lowerx * cos(theta) - blade1lowery * sin(theta),
                                blade1lowerx * sin(theta) + blade1lowery * cos(theta)]
[blade3upperx, blade3uppery] = [blade1upperx * cos(theta) - blade1uppery * sin(theta),
                                blade1upperx * sin(theta) + blade1uppery * cos(theta)]

[x6, y6] = [x0 * cos(-theta) - y0 * sin(-theta), x0 * sin(-theta) + y0 * cos(-theta)]
[x7, y7] = [x1 * cos(-theta) - y1 * sin(-theta), x1 * sin(-theta) + y1 * cos(-theta)]
[x8, y8] = [x2 * cos(-theta) - y2 * sin(-theta), x2 * sin(-theta) + y2 * cos(-theta)]
[x9, y9] = [x3 * cos(-theta) - y3 * sin(-theta), x3 * sin(-theta) + y3 * cos(-theta)]
[x10, y10] = [x4 * cos(-theta) - y4 * sin(-theta), x4 * sin(-theta) + y4 * cos(-theta)]
[x11, y11] = [x5 * cos(-theta) - y5 * sin(-theta), x5 * sin(-theta) + y5 * cos(-theta)]

[x12, y12] = [x6 * cos(-theta) - y6 * sin(-theta), x6 * sin(-theta) + y6 * cos(-theta)]
[x17, y17] = [x11 * cos(-theta) - y11 * sin(-theta), x11 * sin(-theta) + y11 * cos(-theta)]
[x16, y16] = [x10 * cos(-theta) - y10 * sin(-theta), x10 * sin(-theta) + y10 * cos(-theta)]
[x15, y15] = [x9 * cos(-theta) - y9 * sin(-theta), x9 * sin(-theta) + y9 * cos(-theta)]
[x14, y14] = [x8 * cos(-theta) - y8 * sin(-theta), x8 * sin(-theta) + y8 * cos(-theta)]
[x13, y13] = [x7 * cos(-theta) - y7 * sin(-theta), x7 * sin(-theta) + y7 * cos(-theta)]

# point18---p23
xSquare = np.array([-Lsquare / 2, Lsquare / 2, Lsquare / 2, -Lsquare / 2, -Lsquare / 1.6, Lsquare / 1.6])
ySquare = np.array([-Lsquare / 2, -Lsquare / 2, Lsquare / 2, Lsquare / 2, 0, 0])
Gamma = atan((y13 - y6) / (x13 - x6))
[xSquare, ySquare] = [xSquare * cos(Gamma) - ySquare * sin(Gamma), xSquare * sin(Gamma) + ySquare * cos(Gamma)]
[x18, y18] = [xSquare[0], ySquare[0]]
[x19, y19] = [xSquare[1], ySquare[1]]
[x20, y20] = [xSquare[2], ySquare[2]]
[x21, y21] = [xSquare[3], ySquare[3]]
[x22, y22] = [xSquare[4], ySquare[4]]
[x23, y23] = [xSquare[5], ySquare[5]]


[x5a, y5a] = [blade1lowerx[leadregionNumber], blade1lowery[leadregionNumber]]
[x5b, y5b] = [blade1upperx[leadregionNumber], blade1uppery[leadregionNumber]]
[x11a, y11a] = [blade2lowerx[leadregionNumber], blade2lowery[leadregionNumber]]
[x11b, y11b] = [blade2upperx[leadregionNumber], blade2uppery[leadregionNumber]]
[x17a, y17a] = [blade3lowerx[leadregionNumber], blade3lowery[leadregionNumber]]
[x17b, y17b] = [blade3upperx[leadregionNumber], blade3uppery[leadregionNumber]]

x0a = min(
    getInsertPointBetweenCircleAndLine((y0 - y21) / (x0 - x21), -(y0 - y21) * x0 / (x0 - x21) + y0, 0, 0, Rregion))
y0a = (y0 - y21) * (x0a - x0) / (x0 - x21) + y0
[x6a, y6a] = [x0a * cos(-theta) - y0a * sin(-theta), x0a * sin(-theta) + y0a * cos(-theta)]
[x12a, y12a] = [x6a * cos(-theta) - y6a * sin(-theta), x6a * sin(-theta) + y6a * cos(-theta)]

x5c = min(getInsertPointBetweenCircleAndLine(-(x0a - x4) / (y0a - y4),
                                             0.5 * (x0a + x4) * (x0a - x4) / (y0a - y4) + 0.5 * (y0a + y4), 0, 0,
                                             Rregion))
y5c = -(x0a - x4) * x5c / (y0a - y4) + 0.5 * (x0a + x4) * (x0a - x4) / (y0a - y4) + 0.5 * (y0a + y4)
[x11c, y11c] = [x5c * cos(-theta) - y5c * sin(-theta), x5c * sin(-theta) + y5c * cos(-theta)]
[x17c, y17c] = [x11c * cos(-theta) - y11c * sin(-theta), x11c * sin(-theta) + y11c * cos(-theta)]


# caculate 0b
def func(paramlist):
    m, n = paramlist[0], paramlist[1]
    return [(x5 - m) ** 2 + (y5 - n) ** 2 - RendTostart ** 2, (x14 - m) ** 2 + (y14 - n) ** 2 - RendTostart ** 2]


circleCenter = fsolve(func, [0, 0])
# point0b  Perpendicular Bisector
x0b = min(
    getInsertPointBetweenCircleAndLine((y0 - y21) / (x0 - x21), -(y0 - y21) * x0 / (x0 - x21) + y0, circleCenter[0],
                                       circleCenter[1], RendTostart))
y0b = (y0 - y21) * (x0b - x0) / (x0 - x21) + y0
[x6b, y6b] = [x0b * cos(-theta) - y0b * sin(-theta), x0b * sin(-theta) + y0b * cos(-theta)]
[x12b, y12b] = [x6b * cos(-theta) - y6b * sin(-theta), x6b * sin(-theta) + y6b * cos(-theta)]
# test plot
# plt.plot(blade1upperx, blade1uppery, blade1lowerx, blade1lowery)
# plt.plot(blade2upperx, blade2uppery, blade2lowerx, blade2lowery)
# plt.plot(blade3upperx, blade3uppery, blade3lowerx, blade3lowery)
# for i in range(0,24):
#  plt.plot(locals()['x' + str(i)], locals()['y' + str(i)], 'ro', label="point")
# plt.gca().axis('equal')
# plt.plot(x0a,y0a, 'ro', label="point")
# plt.plot(x6a,y6a, 'ro', label="point")
# plt.plot(x12a,y12a, 'ro', label="point")
# plt.plot(x5c,y5c, 'ro', label="point")
# plt.plot(x11c,y11c, 'ro', label="point")
# plt.plot(x17c,y17c, 'ro', label="point")
# plt.plot(x0b,y0b, 'ro', label="point")
# plt.plot(x6b,y6b, 'ro', label="point")
# plt.plot(x12b,y12b, 'ro', label="point")
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
          Vertex(x12, y12, 0, 'v12'),
          Vertex(x13, y13, 0, 'v13'),
          Vertex(x14, y14, 0, 'v14'),
          Vertex(x15, y15, 0, 'v15'),
          Vertex(x16, y16, 0, 'v16'),
          Vertex(x17, y17, 0, 'v17'),
          Vertex(x18, y18, 0, 'v18'),
          Vertex(x19, y19, 0, 'v19'),
          Vertex(x20, y20, 0, 'v20'),
          Vertex(x21, y21, 0, 'v21'),
          Vertex(x22, y22, 0, 'v22'),
          Vertex(x23, y23, 0, 'v23'),
          Vertex(x0a, y0a, 0, 'v0a'),
          Vertex(x6a, y6a, 0, 'v6a'),
          Vertex(x12a, y12a, 0, 'v12a'),
          Vertex(x5a, y5a, 0, 'v5a'),
          Vertex(x5b, y5b, 0, 'v5b'),
          Vertex(x5c, y5c, 0, 'v5c'),
          Vertex(x11a, y11a, 0, 'v11a'),
          Vertex(x11b, y11b, 0, 'v11b'),
          Vertex(x11c, y11c, 0, 'v11c'),
          Vertex(x17a, y17a, 0, 'v17a'),
          Vertex(x17b, y17b, 0, 'v17b'),
          Vertex(x17c, y17c, 0, 'v17c'),
          Vertex(x0b, y0b, 0, 'v0b'),
          Vertex(x6b, y6b, 0, 'v6b'),
          Vertex(x12b, y12b, 0, 'v12b'),
          ]

for v in basevs:
    bmd.add_vertex(v.x, v.y, v.z, v.name + '-z')
    bmd.add_vertex(v.x, v.y, v.z + length_z, v.name + '+z')

# define circle and spline
z = length_z / 2
for i in ['+z', '-z']:
    # outer circle
    midx = max(getInsertPointBetweenCircleAndLine(-(x4 - x3) / (y4 - y3),
                                                  0.5 * (x4 + x3) * (x4 - x3) / (y4 - y3) + 0.5 * (y3 + y4), 0, 0,
                                                  Rregion))
    midy = -(x4 - x3) * midx / (y4 - y3) + 0.5 * (x4 + x3) * (x4 - x3) / (y4 - y3) + 0.5 * (y4 + y3)
    bmd.add_arcedge(('v3' + i, 'v4' + i), 'arc-outer1' + i,
                    Vertex(midx, midy, length_z / 2 + eval(i), 'arc-outer1' + i))
    midx = max(getInsertPointBetweenCircleAndLine(-(x6a - x3) / (y6a - y3),
                                                  0.5 * (x6a + x3) * (x6a - x3) / (y6a - y3) + 0.5 * (y3 + y6a), 0, 0,
                                                  Rregion))
    midy = -(x6a - x3) * midx / (y6a - y3) + 0.5 * (x6a + x3) * (x6a - x3) / (y6a - y3) + 0.5 * (y6a + y3)
    bmd.add_arcedge(('v3' + i, 'v6a' + i), 'arc-outer2' + i,
                    Vertex(midx, midy, length_z / 2 + eval(i), 'arc-outer2' + i))
    midx = max(getInsertPointBetweenCircleAndLine(-(x6a - x11c) / (y6a - y11c),
                                                  0.5 * (x6a + x11c) * (x6a - x11c) / (y6a - y11c) + 0.5 * (y11c + y6a),
                                                  0, 0, Rregion))
    midy = -(x6a - x11c) * midx / (y6a - y11c) + 0.5 * (x6a + x11c) * (x6a - x11c) / (y6a - y11c) + 0.5 * (y6a + y11c)
    bmd.add_arcedge(('v11c' + i, 'v6a' + i), 'arc-outer3' + i,
                    Vertex(midx, midy, length_z / 2 + eval(i), 'arc-outer3' + i))
    midx = max(getInsertPointBetweenCircleAndLine(-(x10 - x11c) / (y10 - y11c),
                                                  0.5 * (x10 + x11c) * (x10 - x11c) / (y10 - y11c) + 0.5 * (y11c + y10),
                                                  0, 0, Rregion))
    midy = -(x10 - x11c) * midx / (y10 - y11c) + 0.5 * (x10 + x11c) * (x10 - x11c) / (y10 - y11c) + 0.5 * (y10 + y11c)
    bmd.add_arcedge(('v11c' + i, 'v10' + i), 'arc-outer4' + i,
                    Vertex(midx, midy, length_z / 2 + eval(i), 'arc-outer4' + i))
    midx = max(getInsertPointBetweenCircleAndLine(-(x10 - x9) / (y10 - y9),
                                                  0.5 * (x10 + x9) * (x10 - x9) / (y10 - y9) + 0.5 * (y9 + y10), 0, 0,
                                                  Rregion))
    midy = -(x10 - x9) * midx / (y10 - y9) + 0.5 * (x10 + x9) * (x10 - x9) / (y10 - y9) + 0.5 * (y10 + y9)
    bmd.add_arcedge(('v9' + i, 'v10' + i), 'arc-outer5' + i,
                    Vertex(midx, midy, length_z / 2 + eval(i), 'arc-outer5' + i))
    midx = min(getInsertPointBetweenCircleAndLine(-(x12a - x9) / (y12a - y9),
                                                  0.5 * (x12a + x9) * (x12a - x9) / (y12a - y9) + 0.5 * (y9 + y12a), 0,
                                                  0, Rregion))
    midy = -(x12a - x9) * midx / (y12a - y9) + 0.5 * (x12a + x9) * (x12a - x9) / (y12a - y9) + 0.5 * (y12a + y9)
    bmd.add_arcedge(('v9' + i, 'v12a' + i), 'arc-outer6' + i,
                    Vertex(midx, midy, length_z / 2 + eval(i), 'arc-outer6' + i))
    midx = min(getInsertPointBetweenCircleAndLine(-(x12a - x17c) / (y12a - y17c),
                                                  0.5 * (x12a + x17c) * (x12a - x17c) / (y12a - y17c) + 0.5 * (
                                                              y17c + y12a), 0, 0, Rregion))
    midy = -(x12a - x17c) * midx / (y12a - y17c) + 0.5 * (x12a + x17c) * (x12a - x17c) / (y12a - y17c) + 0.5 * (
                y12a + y17c)
    bmd.add_arcedge(('v17c' + i, 'v12a' + i), 'arc-outer7' + i,
                    Vertex(midx, midy, length_z / 2 + eval(i), 'arc-outer7' + i))
    midx = min(getInsertPointBetweenCircleAndLine(-(x16 - x17c) / (y16 - y17c),
                                                  0.5 * (x16 + x17c) * (x16 - x17c) / (y16 - y17c) + 0.5 * (y17c + y16),
                                                  0, 0, Rregion))
    midy = -(x16 - x17c) * midx / (y16 - y17c) + 0.5 * (x16 + x17c) * (x16 - x17c) / (y16 - y17c) + 0.5 * (y16 + y17c)
    bmd.add_arcedge(('v17c' + i, 'v16' + i), 'arc-outer8' + i,
                    Vertex(midx, midy, length_z / 2 + eval(i), 'arc-outer8' + i))
    midx = min(getInsertPointBetweenCircleAndLine(-(x16 - x15) / (y16 - y15),
                                                  0.5 * (x16 + x15) * (x16 - x15) / (y16 - y15) + 0.5 * (y15 + y16), 0,
                                                  0, Rregion))
    midy = -(x16 - x15) * midx / (y16 - y15) + 0.5 * (x16 + x15) * (x16 - x15) / (y16 - y15) + 0.5 * (y16 + y15)
    bmd.add_arcedge(('v15' + i, 'v16' + i), 'arc-outer9' + i,
                    Vertex(midx, midy, length_z / 2 + eval(i), 'arc-outer9' + i))
    midx = min(getInsertPointBetweenCircleAndLine(-(x0a - x15) / (y0a - y15),
                                                  0.5 * (x0a + x15) * (x0a - x15) / (y0a - y15) + 0.5 * (y15 + y0a), 0,
                                                  0, Rregion))
    midy = -(x0a - x15) * midx / (y0a - y15) + 0.5 * (x0a + x15) * (x0a - x15) / (y0a - y15) + 0.5 * (y0a + y15)
    bmd.add_arcedge(('v15' + i, 'v0a' + i), 'arc-outer10' + i,
                    Vertex(midx, midy, length_z / 2 + eval(i), 'arc-outer10' + i))
    midx = min(getInsertPointBetweenCircleAndLine(-(x0a - x5c) / (y0a - y5c),
                                                  0.5 * (x0a + x5c) * (x0a - x5c) / (y0a - y5c) + 0.5 * (y5c + y0a), 0,
                                                  0, Rregion))
    midy = -(x0a - x5c) * midx / (y0a - y5c) + 0.5 * (x0a + x5c) * (x0a - x5c) / (y0a - y5c) + 0.5 * (y0a + y5c)
    bmd.add_arcedge(('v5c' + i, 'v0a' + i), 'arc-outer11' + i,
                    Vertex(midx, midy, length_z / 2 + eval(i), 'arc-outer11' + i))
    midx = min(getInsertPointBetweenCircleAndLine(-(x4 - x5c) / (y4 - y5c),
                                                  0.5 * (x4 + x5c) * (x4 - x5c) / (y4 - y5c) + 0.5 * (y5c + y4), 0, 0,
                                                  Rregion))
    midy = -(x4 - x5c) * midx / (y4 - y5c) + 0.5 * (x4 + x5c) * (x4 - x5c) / (y4 - y5c) + 0.5 * (y4 + y5c)
    bmd.add_arcedge(('v5c' + i, 'v4' + i), 'arc-outer12' + i,
                    Vertex(midx, midy, length_z / 2 + eval(i), 'arc-outer12' + i))

    # inner circle
    bmd.add_arcedge(('v0' + i, 'v1' + i), 'arc-inner1' + i,
                    Vertex(Rinner * cos(pi / 2), Rinner * sin(pi / 2), length_z / 2 + eval(i), 'arc-inner1' + i))
    bmd.add_arcedge(('v6' + i, 'v1' + i), 'arc-inner2' + i,
                    Vertex(Rinner * cos(pi / 180), Rinner * sin(pi / 180), length_z / 2 + eval(i), 'arc-inner2' + i))
    bmd.add_arcedge(('v6' + i, 'v7' + i), 'arc-inner3' + i,
                    Vertex(Rinner * cos(-pi / 4), Rinner * sin(-pi / 4), length_z / 2 + eval(i), 'arc-inner3' + i))
    bmd.add_arcedge(('v7' + i, 'v12' + i), 'arc-inner4' + i,
                    Vertex(Rinner * cos(-pi / 2), Rinner * sin(-pi / 2), length_z / 2 + eval(i), 'arc-inner4' + i))
    bmd.add_arcedge(('v12' + i, 'v13' + i), 'arc-inner5' + i,
                    Vertex(Rinner * cos(-5 * pi / 6), Rinner * sin(-5 * pi / 6), length_z / 2 + eval(i),
                           'arc-inner5' + i))
    bmd.add_arcedge(('v13' + i, 'v0' + i), 'arc-inner6' + i,
                    Vertex(Rinner * cos(3 * pi / 4), Rinner * sin(3 * pi / 4), length_z / 2 + eval(i),
                           'arc-inner6' + i))
    # medium circle
    midx = min(getInsertPointBetweenCircleAndLine(-(x14 - x0b) / (y14 - y0b),
                                                  0.5 * (x14 + x0b) * (x14 - x0b) / (y14 - y0b) + 0.5 * (y0b + y14),
                                                  circleCenter[0], circleCenter[1],
                                                  RendTostart))
    midy = -(x14 - x0b) * midx / (y14 - y0b) + 0.5 * (x14 + x0b) * (x14 - x0b) / (y14 - y0b) + 0.5 * (y14 + y0b)
    bmd.add_arcedge(('v0b' + i, 'v14' + i), 'arc-mid-b3-2' + i,
                    Vertex(midx, midy, length_z / 2 + eval(i), 'arc-mid-b3-2' + i))
    bmd.add_arcedge(('v6b' + i, 'v2' + i), 'arc-mid-b1-2' + i,Vertex(midx* cos(-theta) - midy* sin(-theta),midx* sin(-theta) + midy* cos(-theta), length_z / 2 + eval(i), 'arc-mid-b1-2' + i))
    bmd.add_arcedge(('v12b' + i,'v8' + i), 'arc-mid-b2-2' + i,Vertex(midx* cos(theta) - midy* sin(theta),midx* sin(theta) + midy* cos(theta), length_z / 2 + eval(i), 'arc-mid-b2-2' + i))

    midx = min(getInsertPointBetweenCircleAndLine(-(x0b - x5) / (y0b - y5),
                                                  0.5 * (x0b + x5) * (x0b - x5) / (y0b - y5) + 0.5 * (y5 + y0b),
                                                  circleCenter[0], circleCenter[1],
                                                  RendTostart))
    midy = -(x0b - x5) * midx / (y0b - y5) + 0.5 * (x0b + x5) * (x0b - x5) / (y0b - y5) + 0.5 * (y0b + y5)
    bmd.add_arcedge(('v5' + i, 'v0b' + i), 'arc-mid-b1-1' + i,
                    Vertex(midx, midy, length_z / 2 + eval(i), 'arc-mid-b1-1' + i))
    bmd.add_arcedge(('v11' + i, 'v6b' + i), 'arc-mid-b2-1' + i,Vertex(midx* cos(-theta) - midy* sin(-theta),midx* sin(-theta) + midy* cos(-theta), length_z / 2 + eval(i), 'arc-mid-b2-1' + i))
    bmd.add_arcedge(('v17' + i, 'v12b' + i), 'arc-mid-b3-1' + i,Vertex(midx* cos(theta) - midy* sin(theta),midx* sin(theta) + midy* cos(theta), length_z / 2 + eval(i), 'arc-mid-b3-1' + i))


# blades
    upperpoints1 = []
    for j in range(1,leadregionNumber):
        upperpoints1.append(Point(blade1upperx[j],blade1uppery[j], length_z / 2 + eval(i)))
    lowerpoints1 = []
    for j in range(1,leadregionNumber):
        lowerpoints1.append(Point(blade1lowerx[j],blade1lowery[j], length_z / 2 + eval(i)))
    bmd.add_splineedge(('v5' + i, 'v5b' + i), 'blade1-p1' + i, upperpoints1)
    bmd.add_splineedge(('v5' + i, 'v5a' + i), 'blade1-m1' + i, lowerpoints1)

    upperpoints2 = []
    for j in range(leadregionNumber+ 1, len(x)):
        upperpoints2.append(Point(blade1upperx[j], blade1uppery[j], length_z / 2 + eval(i)))
    lowerpoints2 = []
    for j in range(leadregionNumber+ 1, len(x)):
        lowerpoints2.append(Point(blade1lowerx[j], blade1lowery[j], length_z / 2 + eval(i)))
    bmd.add_splineedge(('v5b' + i, 'v2' + i), 'blade1-p2' + i, upperpoints2)
    bmd.add_splineedge(('v5a' + i, 'v2' + i), 'blade1-m2' + i, lowerpoints2)


    upperpoints1 = []
    for j in range(1,leadregionNumber):
        upperpoints1.append(Point(blade2upperx[j],blade2uppery[j], length_z / 2 + eval(i)))
    lowerpoints1 = []
    for j in range(1,leadregionNumber):
        lowerpoints1.append(Point(blade2lowerx[j],blade2lowery[j], length_z / 2 + eval(i)))
    bmd.add_splineedge(('v11' + i, 'v11b' + i), 'blade2-p1' + i, upperpoints1)
    bmd.add_splineedge(('v11' + i, 'v11a' + i), 'blade2-m1' + i, lowerpoints1)

    upperpoints2 = []
    for j in range(leadregionNumber+ 1, len(x)):
        upperpoints2.append(Point(blade2upperx[j], blade2uppery[j], length_z / 2 + eval(i)))
    lowerpoints2 = []
    for j in range(leadregionNumber+ 1, len(x)):
        lowerpoints2.append(Point(blade2lowerx[j], blade2lowery[j], length_z / 2 + eval(i)))
    bmd.add_splineedge(('v11b' + i, 'v8' + i), 'blade2-p2' + i, upperpoints2)
    bmd.add_splineedge(('v11a' + i, 'v8' + i), 'blade2-m2' + i, lowerpoints2)



    upperpoints1 = []
    for j in range(1,leadregionNumber):
        upperpoints1.append(Point(blade3upperx[j],blade3uppery[j], length_z / 2 + eval(i)))
    lowerpoints1 = []
    for j in range(1,leadregionNumber):
        lowerpoints1.append(Point(blade3lowerx[j],blade3lowery[j], length_z / 2 + eval(i)))
    bmd.add_splineedge(('v17' + i, 'v17b' + i), 'blade3-p1' + i, upperpoints1)
    bmd.add_splineedge(('v17' + i, 'v17a' + i), 'blade3-m1' + i, lowerpoints1)

    upperpoints2 = []
    for j in range(leadregionNumber+ 1, len(x)):
        upperpoints2.append(Point(blade3upperx[j], blade3uppery[j], length_z / 2 + eval(i)))
    lowerpoints2 = []
    for j in range(leadregionNumber+ 1, len(x)):
        lowerpoints2.append(Point(blade3lowerx[j], blade3lowery[j], length_z / 2 + eval(i)))
    bmd.add_splineedge(('v17b' + i, 'v14' + i), 'blade3-p2' + i, upperpoints2)
    bmd.add_splineedge(('v17a' + i, 'v14' + i), 'blade3-m2' + i, lowerpoints2)


#
def vnamegen(x0z0, x1z0, x1z1, x0z1):
    return (x0z0 + '-z', x1z0 + '-z', x1z1 + '-z', x0z1 + '-z',
            x0z0 + '+z', x1z0 + '+z', x1z1 + '+z', x0z1 + '+z')



# == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
# create blocks
# 123456789: to be replaced
b0 = bmd.add_hexblock(vnamegen('v0', 'v1', 'v2', 'v5a'), (FCells, DCells, 1), 'b0',
                      SimpleGrading(123456789, 1 / DGrading, 1))
b1 = bmd.add_hexblock(vnamegen('v5b', 'v2', 'v3', 'v4'), (FCells, UpCells, 1), 'b1',
                      SimpleGrading(123456789, UpGrading, 1))
b2 = bmd.add_hexblock(vnamegen('v1', 'v6', 'v6b', 'v2'), (FCells, DCells, 1), 'b2',
                      SimpleGrading(123456789, 1 / DGrading, 1))
b3 = bmd.add_hexblock(vnamegen('v2', 'v6b', 'v6a', 'v3'), (FCells, UpCells, 1), 'b3',
                      SimpleGrading(123456789, UpGrading, 1))
b4 = bmd.add_hexblock(vnamegen('v6', 'v7', 'v8', 'v11a'), (FCells, DCells, 1), 'b4',
                      SimpleGrading(123456789, 1 / DGrading, 1))
b5 = bmd.add_hexblock(vnamegen('v11b', 'v8', 'v9', 'v10'), (FCells, UpCells, 1), 'b5',
                      SimpleGrading(123456789, UpGrading, 1))
b6 = bmd.add_hexblock(vnamegen('v7', 'v12', 'v12b', 'v8'), (FCells, DCells, 1), 'b6',
                      SimpleGrading(123456789, 1 / DGrading, 1))
b7 = bmd.add_hexblock(vnamegen('v8', 'v12b', 'v12a', 'v9'), (FCells, UpCells, 1), 'b7',
                      SimpleGrading(123456789, UpGrading, 1))
b8 = bmd.add_hexblock(vnamegen('v12', 'v13', 'v14', 'v17a'), (FCells, DCells, 1), 'b8',
                      SimpleGrading(123456789, 1 / DGrading, 1))
b9 = bmd.add_hexblock(vnamegen('v17b', 'v14', 'v15', 'v16'), (FCells, UpCells, 1), 'b9',
                      SimpleGrading(123456789, UpGrading, 1))
b10 = bmd.add_hexblock(vnamegen('v13', 'v0', 'v0b', 'v14'), (FCells, DCells, 1), 'b10',
                       SimpleGrading(123456789, 1 / DGrading, 1))
b11 = bmd.add_hexblock(vnamegen('v14', 'v0b', 'v0a', 'v15'), (FCells, UpCells, 1), 'b11',
                       SimpleGrading(123456789, UpGrading, 1))
b12 = bmd.add_hexblock(vnamegen('v21', 'v20', 'v1', 'v0'), (FCells, CCells, 1), 'b12',
                       SimpleGrading(123456789, 1 / CGrading, 1))
b13 = bmd.add_hexblock(vnamegen('v20', 'v23', 'v6', 'v1'), (FCells, CCells, 1), 'b13',
                       SimpleGrading(123456789, 1 / CGrading, 1))
b14 = bmd.add_hexblock(vnamegen('v23', 'v19', 'v7', 'v6'), (FCells, CCells, 1), 'b14',
                       SimpleGrading(123456789, 1 / CGrading, 1))
b15 = bmd.add_hexblock(vnamegen('v19', 'v18', 'v12', 'v7'), (FCells, CCells, 1), 'b15',
                       SimpleGrading(123456789, 1 / CGrading, 1))
b16 = bmd.add_hexblock(vnamegen('v18', 'v22', 'v13', 'v12'), (FCells, CCells, 1), 'b16',
                       SimpleGrading(123456789, 1 / CGrading, 1))
b17 = bmd.add_hexblock(vnamegen('v22', 'v21', 'v0', 'v13'), (FCells, CCells, 1), 'b17',
                       SimpleGrading(123456789, 1 / CGrading, 1))
b18 = bmd.add_hexblock(vnamegen('v22', 'v23', 'v20', 'v21'), (FCells, FCells, 1), 'b18',
                       SimpleGrading(123456789, 123456789, 1))
b19 = bmd.add_hexblock(vnamegen('v23', 'v22', 'v18', 'v19'), (FCells, FCells, 1), 'b19',
                       SimpleGrading(123456789, 123456789, 1))
b5a = bmd.add_hexblock(vnamegen('v0', 'v5a', 'v5', 'v0b'), (DCells, DCells, 1), 'b5a',
                       SimpleGrading(1 / DGrading, 1 / DGrading, 1))
b5b = bmd.add_hexblock(vnamegen('v5', 'v5b', 'v4', 'v5c'), (DCells, UpCells, 1), 'b5b',
                       SimpleGrading(DGrading, UpGrading, 1))
b5c = bmd.add_hexblock(vnamegen('v0b', 'v5', 'v5c', 'v0a'), (DCells, UpCells, 1), 'b5c',
                       SimpleGrading(1 / DGrading, UpGrading, 1))
b11a = bmd.add_hexblock(vnamegen('v6', 'v11a', 'v11', 'v6b'), (DCells, DCells, 1), 'b11a',
                        SimpleGrading(1 / DGrading, 1 / DGrading, 1))
b11b = bmd.add_hexblock(vnamegen('v11', 'v11b', 'v10', 'v11c'), (DCells, UpCells, 1), 'b11b',
                        SimpleGrading(DGrading, UpGrading, 1))
b11c = bmd.add_hexblock(vnamegen('v6b', 'v11', 'v11c', 'v6a'), (DCells, UpCells, 1), 'b11c',
                        SimpleGrading(1 / DGrading, UpGrading, 1))
b17a = bmd.add_hexblock(vnamegen('v12', 'v17a', 'v17', 'v12b'), (DCells, DCells, 1), 'b17a',
                        SimpleGrading(1 / DGrading, 1 / DGrading, 1))
b17b = bmd.add_hexblock(vnamegen('v17', 'v17b', 'v16', 'v17c'), (DCells, UpCells, 1), 'b17b',
                        SimpleGrading(DGrading, UpGrading, 1))
b17c = bmd.add_hexblock(vnamegen('v12b', 'v17', 'v17c', 'v12a'), (DCells, UpCells, 1), 'b17c',
                        SimpleGrading(1 / DGrading, UpGrading, 1))
# == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
# boundary
bmd.add_boundary('wall', 'blade1-l2', [b0.face('yp')])
bmd.add_boundary('wall', 'blade1-h2', [b1.face('ym')])
bmd.add_boundary('wall', 'blade1-l1', [b5a.face('xp')])
bmd.add_boundary('wall', 'blade1-h1', [b5b.face('ym')])

bmd.add_boundary('wall', 'blade2-l2', [b4.face('yp')])
bmd.add_boundary('wall', 'blade2-h2', [b5.face('ym')])
bmd.add_boundary('wall', 'blade2-l1', [b11a.face('xp')])
bmd.add_boundary('wall', 'blade2-h1', [b11b.face('ym')])

bmd.add_boundary('wall', 'blade3-l2', [b8.face('yp')])
bmd.add_boundary('wall', 'blade3-h2', [b9.face('ym')])
bmd.add_boundary('wall', 'blade3-l1', [b17a.face('xp')])
bmd.add_boundary('wall', 'blade3-h1', [b17b.face('ym')])

bmd.add_boundary('patch', 'AMIin-blade1-1', [b5c.face('yp')])
bmd.add_boundary('patch', 'AMIin-blade1-2', [b5b.face('yp')])
bmd.add_boundary('patch', 'AMIin-blade1-3', [b1.face('yp')])
bmd.add_boundary('patch', 'AMIin-blade1-4', [b3.face('yp')])

bmd.add_boundary('patch', 'AMIin-blade2-1', [b11c.face('yp')])
bmd.add_boundary('patch', 'AMIin-blade2-2', [b11b.face('yp')])
bmd.add_boundary('patch', 'AMIin-blade2-3', [b5.face('yp')])
bmd.add_boundary('patch', 'AMIin-blade2-4', [b7.face('yp')])

bmd.add_boundary('patch', 'AMIin-blade3-1', [b17c.face('yp')])
bmd.add_boundary('patch', 'AMIin-blade3-2', [b17b.face('yp')])
bmd.add_boundary('patch', 'AMIin-blade3-3', [b9.face('yp')])
bmd.add_boundary('patch', 'AMIin-blade3-4', [b11.face('yp')])



# == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
# prepare for output
bmd.assign_vertexid()
f = open('blockMeshDict', 'wb')
f.write(bmd.format().encode())
f.close()
with fileinput.FileInput('blockMeshDict', inplace=True) as file:
    for line in file:
        print(line.replace(str(123456789), FGrading), end='')
