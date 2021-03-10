import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
import matplotlib.colors as mcolors
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import os
import re
import numpy as np
import h5py
import sys
from os.path import dirname, realpath
pypath = dirname(dirname(dirname(realpath(__file__)))) + '/python/'
sys.path.insert(1, pypath)
import utils
from odes import car2

dirpath = dirname(realpath(__file__))

transfile = '/abstfull_0.2-0.2-0.2.h5'
labelfile = '/labels_dba1_abstfull_0.2-0.2-0.2.h5'
ctlrfile = '/controller_dba1_0.2-0.2-0.2.h5'
specfile = '/dba1.txt'


# # Simulation of static motion planning
tau, X, eta, _, winids, controller = \
    utils.read_controller_abst_from_h5(dirpath+ctlrfile)
# # Compute the percentage of winning set on the state space
winset = controller.xgrid[winids, :]
print("\nWinning set coverage:")
winper = "{:.2%}".format(winids.size/controller.xgrid.shape[0])
print(winper)
# # Load specification
dba = utils.read_spec_from_txt(dirpath+specfile)
# # Simulation
Tsim = 50
num_acc = 3
x0 = np.array([1.0, 1.0, np.pi/3.])
i0 = utils.index_in_grid(x0, controller.xgrid)
if(not np.any(winids == i0)):
    sys.exit("The initial condition is not in the winning set.\n")
xsim, usim, qsim, tsim = utils.simulate_abstbased_dba_control(
    tau, Tsim, num_acc, x0, car2, dba, controller)


# # Display workspace
fig = plt.figure()
ax = plt.axes()
ax.set_xlim(0,10)
ax.set_ylim(0,10)

xgrid = np.array([])
eta = np.array([])
labels = np.array([])
obs = np.array([])
with h5py.File(dirpath+transfile, 'r') as ft,\
     h5py.File(dirpath+labelfile, 'r') as fl:
    eta = ft['eta'][...]
    xgrid = ft['xgrid'][...]
    obs = ft['obs'][...]
    labels = fl['labels'][...]

oset = xgrid[obs, :]
ax.add_collection(
    utils.polycoll_grid_array(oset, eta, True, 'gray', 0.7)
)
gset = xgrid[np.where(labels>0),:].squeeze()
ax.add_collection(
    utils.polycoll_grid_array(gset, eta, True, 'palegreen', 0.7)
)

obstacles = [[0.0, 0.5, 5.0, 6.0],
             [2.4, 2.6, 0.0, 3.2],
             [3.9, 4.1, 9.0, 10.0], #[3.9, 4.1, 8.0, 10.0]
             [5.9, 6.1, 0.0, 0.6],
             [5.9, 6.1, 3.8, 5.1], # [5.9, 6.1, 3.8, 6.1],
             [6.1, 10.0, 4.9, 5.1]] # [6.1, 10.0, 5.9, 6.1]
rects_obs = [patches.Rectangle((obstacles[i][0], obstacles[i][2]),
                               obstacles[i][1]-obstacles[i][0],
                               obstacles[i][3]-obstacles[i][2],
                               linewidth=1,edgecolor='k',facecolor='k')
              for i in range(len(obstacles))]

goals = [[0.5, 2.0, 7.5, 9.5],  # a
         [7.5, 9.5, 0.8, 3.0]]  # d
rects_gs = [patches.Rectangle((goals[i][0], goals[i][2]),
                               goals[i][1]-goals[i][0],
                               goals[i][3]-goals[i][2],
                               linewidth=1,edgecolor='y',fill=False)
            for i in range(len(goals))]
circ_gs = patches.Circle([8.0, 8.0], radius=0.8, linewidth=1,
                         edgecolor='y', fill=False)


for rect in rects_obs:
    ax.add_patch(rect)

for rect in rects_gs:
    ax.add_patch(rect)
ax.add_patch(circ_gs)

ax.text((goals[0][0]+goals[0][1])/2.0-0.5, (goals[0][2]+goals[0][3])/2.0, "pickup")
ax.text((goals[1][0]+goals[1][1])/2.0-0.5, (goals[1][2]+goals[1][3])/2.0, "drop")
ax.text(7.8, 7.8, "count")


# # # Plot winning set in 2D plane
# ax.add_collection(
#     utils.polycoll_grid_array(winset, eta, True, 'palegoldenrod', 0.7)
# )


# # Plot closed-loop trajectoryies
colors = list(mcolors.TABLEAU_COLORS.keys())
qdiff = np.diff(np.insert(qsim, 0, dba.q0))
qchks = np.argwhere(qdiff != 0).squeeze()
qchks = np.insert(qchks, 0, 0)
qchks = np.append(qchks, qsim.size)
q = dba.q0
for k in range(qchks.size-1):
    q = q + qdiff[qchks[k]]
    xs = xsim[qchks[k]:qchks[k+1]+1, :]
    ax.plot(xs[:, 0], xs[:, 1], color=colors[q])

ax.plot(xsim[0, 0], xsim[0, 1], marker='^',
        markerfacecolor='r', markeredgecolor='r')
ax.plot(xsim[-1, 0], xsim[-1, 1], marker='v',
        markerfacecolor='g', markeredgecolor='g')



plt.show()
