# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Add the SPLINTER directory to the search path, so we can include it
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from os import sys, path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import splinterpy

# Only for dev purposes
import os.path
if os.path.isdir("/home/bjarne/"):
    splinterpy.load("/home/bjarne/Code/splinter/build/debug/libsplinter-4-0.so")
elif os.path.isdir("/home/anders/"):
    splinterpy.load("/home/anders/SPLINTER/build/debug/libsplinter-4-0.so")


# Example with two variables
def f(x):
    return x[0]*x[1]


x1 = np.arange(-5, 5, 0.25)
x2 = np.arange(-5, 5, 0.25)
X1, X2 = np.meshgrid(x1, x2)
Y = np.sqrt(X1 ** 2 + X2 ** 2)

data = []
x = []
y = []
for i in range(len(x1)):
    for j in range(len(x2)):
        x1_ij = X1[j, i]
        x2_ij = X2[j, i]

        x.append([x1_ij, x2_ij])
        y.append(Y[i, j])

# Cubic B-spline
bspline = splinterpy.bspline_interpolator(x, y, degree=3)

Zbs = Y

for i in range(len(x1)):
    for j in range(len(x2)):
        x1_ij = X1[i, j]
        x2_ij = X2[i, j]
        Zbs[i, j] = bspline.eval([x1_ij, x2_ij])

# Plot f
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X1, X2, Y, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
# ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
fig.colorbar(surf, shrink=0.5, aspect=5)

# Plot b-spline
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X1, X2, Zbs, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
# ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
fig.colorbar(surf, shrink=0.5, aspect=5)

# TODO: plot grid using plot_wireframe()

plt.show()
