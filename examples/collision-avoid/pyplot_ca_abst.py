import h5py
import numpy as np
import re
from functools import reduce
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PolyCollection
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from os.path import dirname, realpath
import argparse
parser = argparse.ArgumentParser(description='Indicate control scenarios.')
parser.add_argument('suffix', type=str, help='Indicate the suffix of the file')
args = parser.parse_args()


dirpath = dirname(realpath(__file__))

filename = '/controller_safety_abst-' + args.suffix + '.h5'
rmin = 0.8

f = h5py.File(dirpath+filename, "r")
tau = f['ts'][...]
goalset = f['G'][...]
statespace = f['X'][...]
inputspace = f['U'][...]
xgrid = f['xgrid'][...]
winids = f['WinSet'][...]
iswin = f['Win'][...]
ctlr = f['LeastCtlr'][...]
f.close()

tau = tau[0]
eta = np.array([0.2, 0.2, 0.2])
vmin = xgrid[0,:] - eta/2.
vmax = xgrid[-1,:] + eta/2.
numid_dim = (xgrid[-1,:] - xgrid[0,:])/eta + 1
numid_dim = numid_dim.astype(int)

def val_to_id(x, numid_dim, vmin, eta):
    ids = np.round((x - vmin)/eta).astype(int)
    return (ids[2]*numid_dim[1]+ids[1])*numid_dim[0]+ids[0]

def is_inside(x, w):
    return np.all(x>w[:,0]) and np.all(x<w[:,1])

# # 3D plot of center of winning grids
# w_center = np.zeros((winset.shape[0], 3))
# w_center[:, 0] = (winset[:, 0]+winset[:, 1])/2.0
# w_center[:, 1] = (winset[:, 2]+winset[:, 3])/2.0
# w_center[:, 2] = (winset[:, 4]+winset[:, 5])/2.0

# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.scatter(w_center[:, 0], w_center[:, 1], w_center[:, 2])

# # x-y 2D plot of winning set
winset = xgrid[winids[:-1], :]
xys = winset[:, 0:2]
tmp = np.zeros((xys.shape[0],8))
tmp[:, 0] = xys[:, 0] - eta[0]/2.
tmp[:, 1] = xys[:, 1] - eta[1]/2.
tmp[:, 2] = xys[:, 0] + eta[0]/2.
tmp[:, 3] = xys[:, 1] - eta[1]/2.
tmp[:, 4] = xys[:, 0] + eta[0]/2.
tmp[:, 5] = xys[:, 1] + eta[1]/2.
tmp[:, 6] = xys[:, 0] - eta[0]/2.
tmp[:, 7] = xys[:, 1] + eta[1]/2.
verts = tmp.reshape(xys.shape[0], 4, 2)

fig, ax = plt.subplots()
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)

coll = PolyCollection(verts,closed=True,edgecolors='k',facecolors='g',alpha=0.5)
ax.add_collection(coll)

circ = patches.Circle((0., 0.), radius=rmin, fill=False, color='r')
ax.add_patch(circ)

# # Simulation
def twoagent(t, x, u, d):
    return np.array([-u[0]+d[0]*np.cos(x[2])+u[1]*x[1], d[0]*np.sin(x[2])-u[1]*x[0], d[1]-u[1]])

rng = np.random.default_rng()

Tsim = 50
x0 = np.array([-1.3, 1.3, 0])

t = 0
x = x0
tsim = []
xsim = []
usim = []
dsim = []
while(t<Tsim and is_inside(x, statespace)):
    i = val_to_id(x, numid_dim, vmin, eta) #convert x to node id
    if(iswin[i]):
        uset = np.argwhere(ctlr[i, :]).squeeze() #get all valid inputs
        if(len(uset) > 1):
            uid = rng.choice(uset, 1) #randomly pick one
        else:
            uid = uset[0]
        u = inputspace[uid, :].squeeze()
        d = np.random.uniform(-0.8, 0.8, 2) #generate random disturbace
        # Integrate ode
        sol = solve_ivp(twoagent, [0, tau], x, method='RK45', args=(u, d))
        tt = sol.t[-1]
        y = sol.y[:,-1]
        if(y[2] > np.pi):
            y[2] -= 2*np.pi
        if(y[2] < -np.pi):
            y[2] += 2*np.pi
        # Save trajectories
        tsim.append(t)
        xsim.append(x)
        usim.append(u)
        dsim.append(d)
        # Update state
        x = y
        t += tt
    else:
        print("State is not inside the winning set.\n")
        break

xsim = np.asarray(xsim)
ax.plot(xsim[:, 0], xsim[:, 1], 'k')
ax.plot(xsim[0, 0], xsim[0, 1], marker='o', markerfacecolor='r')


# # # Load trajectories
# egofile = "/Users/yinan/Desktop/rocs/examples/collision-avoid/sim_closedloop_ca.txt"
# otherfile = "/Users/yinan/Desktop/rocs/examples/collision-avoid/sim_other_robot.txt"
# x = []
# u = []
# xo = []
# with open(egofile, "r") as ego, open(otherfile, "r") as other:
#     lines = ego.readlines()
#     for line in lines:
#         line = line.strip()
#         xu = re.split(';', line)
#         x.append([float(e) for e in re.split(',', xu[0])])
#         u.append([float(e) for e in re.split(',|\n', xu[1])])
#     olines = other.readlines()
#     for line in olines:
#         line = line.strip()
#         xo.append([float(e) for e in re.split(',', line)])

# x = np.array(x)
# u = np.array(u)
# xo = np.array(xo)
# xr = xo - x
# ax.plot(xr[76:-1, 0], xr[76:-1, 1], 'k')
# ax.plot(xr[76, 0], xr[76, 1], marker='v',markerfacecolor='b')
# ax.plot(xr[-1, 0], xr[-1, 1], marker='^',markerfacecolor='m')

plt.show()
