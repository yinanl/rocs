import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import h5py

import numpy as np
import re

filename = "/Users/yinan/Desktop/rocs/examples/collision-avoid/\
controller_safety_itvl_0.8-1.2-0.3-0.2.h5"

f = h5py.File(filename, "r")
keys = list(f.keys())
goalset = f['G'][...]
tag = f['tag'][...]
pavings = f['pavings'][...]
winset = pavings[np.argwhere(tag == 1).squeeze()]
f.close()

# # 3D plot of center of winning grids
# w_center = np.zeros((winset.shape[0], 3))
# w_center[:, 0] = (winset[:, 0]+winset[:, 1])/2.0
# w_center[:, 1] = (winset[:, 2]+winset[:, 3])/2.0
# w_center[:, 2] = (winset[:, 4]+winset[:, 5])/2.0

# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.scatter(w_center[:, 0], w_center[:, 1], w_center[:, 2])

# # x-y 2D plot of winning set
xys = winset[:, 0:4]
tmp = np.zeros((xys.shape[0],xys.shape[1]*2))
tmp[:, 0] = xys[:, 0]
tmp[:, 1] = xys[:, 2]
tmp[:, 2] = xys[:, 1]
tmp[:, 3] = xys[:, 2]
tmp[:, 4] = xys[:, 1]
tmp[:, 5] = xys[:, 3]
tmp[:, 6] = xys[:, 0]
tmp[:, 7] = xys[:, 3]
verts = tmp.reshape(xys.shape[0], 4, 2)

fig, ax = plt.subplots()
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)

coll = PolyCollection(verts,closed=True,edgecolors='k',facecolors='g',alpha=0.5)
ax.add_collection(coll)

# # Load trajectories
egofile = "/Users/yinan/Desktop/rocs/examples/collision-avoid/sim_closedloop_ca.txt"
otherfile = "/Users/yinan/Desktop/rocs/examples/collision-avoid/sim_other_robot.txt"
x = []
u = []
xo = []
with open(egofile, "r") as ego, open(otherfile, "r") as other:
    lines = ego.readlines()
    for line in lines:
        line = line.strip()
        xu = re.split(';', line)
        x.append([float(e) for e in re.split(',', xu[0])])
        u.append([float(e) for e in re.split(',|\n', xu[1])])
    olines = other.readlines()
    for line in olines:
        line = line.strip()
        xo.append([float(e) for e in re.split(',', line)])

x = np.array(x)
u = np.array(u)
xo = np.array(xo)
xr = xo - x
ax.plot(xr[76:-1, 0], xr[76:-1, 1], 'k')
ax.plot(xr[76, 0], xr[76, 1], marker='v',markerfacecolor='b')
ax.plot(xr[-1, 0], xr[-1, 1], marker='^',markerfacecolor='m')

plt.show()
