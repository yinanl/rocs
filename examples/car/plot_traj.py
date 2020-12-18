import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
import os
import re
import numpy as np


# # load simulated trajectories
folder = "/Users/yinan/Desktop/rocs/examples/car/"
spec = "dba3"

filename = folder + "sim_traj_" + spec + ".txt"
x = []
u = []
xo = []
with open(filename, "r") as ego:
    lines = ego.readlines()
    for line in lines:
        line = line.strip()
        xu = re.split(';', line)
        x.append([float(e) for e in re.split(',', xu[0])])
        u.append([float(e) for e in re.split(',|\n', xu[1])])

x = np.array(x)
u = np.array(u)


# # Set up workspace
theta = 3.5
nG = 3
nA = 4
G = np.zeros(shape=(3, 2, 3))
G[:, :, 0] = np.array([[1.0, 2.0], [0.5, 2.0], [-theta, theta]])
G[:, :, 1] = np.array([[0.5, 2.5], [7.5, 8.5], [-theta, theta]])
G[:, :, 2] = np.array([[7.1, 9.1], [4.6, 6.4], [-theta, theta]])
A = np.zeros(shape=(3, 2, 4))
A[:, :, 0] = np.array([[1.6, 5.7], [4.0, 5.0], [-theta, theta]])
A[:, :, 1] = np.array([[3.0, 5.0], [5.0, 8.0], [-theta, theta]])
A[:, :, 2] = np.array([[4.3, 5.7], [1.8, 4.0], [-theta, theta]])
A[:, :, 3] = np.array([[5.7, 8.5], [1.8, 2.5], [-theta, theta]])

X = np.array([[0, 10], [0, 10], [-theta, theta]])

# # x-y 2D plot of winning set
fig, ax = plt.subplots()
ax.set_xlim(X[0, 0], X[0, 1])
ax.set_ylim(X[1, 0], X[1, 1])
FS = 16
LW = 2

for i in range(nG):

    ax.add_patch(
        patches.Rectangle((G[0, 0, i], G[1, 0, i]),
                          G[0, 1, i]-G[0, 0, i], G[1, 1, i]-G[1, 0, i],
                          facecolor='gold')
    )

for i in range(nA):
    ax.add_patch(
        patches.Rectangle((A[0, 0, i], A[1, 0, i]),
                          A[0, 1, i]-A[0, 0, i], A[1, 1, i]-A[1, 0, i],
                          facecolor='dimgray')
    )
plt.text(1.4, 1.2, '$a_1$', fontsize=FS)
plt.text(1.3, 8.0, '$a_2$', fontsize=FS)
plt.text(8.0, 5.5, '$a_3$', fontsize=FS)


# # Plot closed-loop trajectory
ax.plot(x[:, 0], x[:, 1])
ax.plot(x[0,0], x[0,1], marker='^', markerfacecolor='r')
ax.plot(x[-1,0], x[-1,1], marker='v', markerfacecolor='g')


plt.show()
