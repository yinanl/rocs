import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
from scipy.integrate import solve_ivp
import re
import numpy as np
import math
import h5py
import argparse
import sys
from os.path import dirname, realpath
pypath = dirname(dirname(dirname(realpath(__file__)))) + '/python/'
sys.path.insert(1, pypath)
import utils
from odes import car2


dirpath = dirname(realpath(__file__))


# # # Test local reachavoid
# node = []
# with open(dirpath+'logs_subgraph.txt', "r") as sub:
#     lines = sub.readlines()
#     for line in lines:
#         line = line.strip()
#         node.append(int(line))
# node = np.array(node)

# # print(np.intersect1d(avoid, node))
# unique = np.unique(node)
# unique.shape


# # Display workspace
goalclr = 'lightgreen'
gclr = 'olivedrab'
lfont = {'fontname':'Consolas'}
tfont = {'fontname':'Times New Roman'}

fig = plt.figure()
ax = plt.axes()
ax.set_xlim(0,10)
ax.set_ylim(0,10)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_title('2D Collision-free Trajectory', **tfont)


transfile = '/abstfull_0.2-0.2-0.2.h5'
ptfile = '/logs_replan.h5'
xgrid = np.array([])
eta = np.array([])
obs = np.array([])
avoid = np.array([])
target = np.array([])
x = np.array([])
xo = np.array([])
with h5py.File(dirpath+transfile, 'r') as ft,\
     h5py.File(dirpath+ptfile, 'r') as f:
    eta = ft['eta'][...]
    xgrid = ft['xgrid'][...]
    obs = ft['obs'][...]
    avoid = f['avoid'][...]
    target = f['target'][...]
    x = f['x'][...]
    xo = f['xo'][...]

oset = xgrid[obs, :]
ax.add_collection(
    utils.polycoll_grid_array(oset, eta, True, 'k', 0.7)
)
aset = xgrid[avoid, :]
ax.add_collection(
    utils.polycoll_grid_array(aset, eta, True, 'gray', 0.7)
)
gset = xgrid[target, :]
ax.add_collection(
    utils.polycoll_grid_array(gset, eta, True, goalclr, 0.7)
)

# xo = np.array([3.08789, 5.67336, -1.31934])
ax.plot(xo[0], xo[1], marker='o', markeredgecolor='r',markerfacecolor='r')
# x0 = np.array([4.02228, 5.45505, 3.0572])
ax.plot(x[0], x[1], marker='o', markeredgecolor='b', markerfacecolor='b')
# # Two target centers
zfile = '/logs_newgoals.txt'
z = []
with open(dirpath+zfile, "r") as f:
    lines = f.readlines()
    for line in lines:
        line = line.strip()
        z.append([float(e) for e in re.split(' ', line)])
z = np.array(z)
z1 = z[0, :]
z2 = z[1, :]
ax.plot(z1[0], z1[1], marker='o', markeredgecolor='olivedrab',
        markerfacecolor='olivedrab')
ax.plot(z2[0], z2[1], marker='o', markeredgecolor='olivedrab',
        markerfacecolor='olivedrab')


# # Simulate local reach control
gfile = '/controller_dba1_0.2-0.2-0.2.h5'
tau, X, eta, _, winids, controller = \
    utils.read_controller_abst_from_h5(dirpath+gfile)

specfile = '/dba1.txt'
dba = utils.read_spec_from_txt(dirpath+specfile)

ctlrfile = 'logs_newctlr.txt'
ctlr = []
with open(dirpath+ctlrfile, "r") as f:
    lines = f.readlines()
    for line in lines:
        line = line.strip()
        ctlr.append([float(e) for e in re.split(',', line)])
ctlr = np.array(ctlr)

rng = np.random.default_rng()
q = 0
t = 0
x = x0
tsim = []
xsim = []
usim = []
qsim = []
while(True):
    i = utils.index_in_grid(x, xgrid)  # convert x to node id
    xp = i*dba.n_dba + q # (n0xn1)-based product id
    uid = ctlr[np.argwhere(ctlr[:,0]==xp), 1].astype(int)
    if(uid.size==0):
        break;
    u = controller.ugrid[uid, :].squeeze()

    # Integrate ode
    sol = solve_ivp(car2, [0, tau], x, method='RK45', args=(u,))
    tt = sol.t[-1]
    y = sol.y[:, -1]
    if(y[2] > np.pi):
        y[2] -= 2*np.pi
    if(y[2] < -np.pi):
        y[2] += 2*np.pi

    # Save trajectories
    tsim.append(t)
    xsim.append(x)
    usim.append(u)
    qsim.append(q)

    # Update state
    x = y
    t += tt

xsim = np.asarray(xsim)
tsim = np.asarray(tsim)
usim = np.asarray(usim)
qsim = np.asarray(qsim)

ax.plot(xsim[:, 0], xsim[:, 1])

plt.show()
