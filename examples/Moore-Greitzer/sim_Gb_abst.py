import h5py
import numpy as np
import os

from scipy.integrate import odeint
from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PolyCollection
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection


dir_path = os.path.dirname(os.path.realpath(__file__))
filename = "controller_Gb_abst_0.00018.h5"

f = h5py.File(dir_path+filename, "r")
tau = f['ts'][...]
goalset = f['G'][...]
X = f['X'][...]
U = f['U'][...]
xgrid = f['xgrid'][...]
winids = f['WinSet'][...]
eta = f['eta'][...]
ctlr = f['OptCtlr'][...]
encode3 = f['encode3'][...]
nts_ctrlr = f['nts_ctrlr'][...]
q_prime = f['q_prime'][...]
f.close()

tau = tau[0]
goal_real = np.array([[0.5009, 0.5069], [0.6575, 0.6635]])

nDBA = 2
q0 = 0


vmin = xgrid[0,:] - eta/2.
vmax = xgrid[-1,:] + eta/2.
numid_dim = (xgrid[-1,:] - xgrid[0,:])/eta + 1
numid_dim = numid_dim.astype(int)


def val_to_id(x, numid_dim, vmin, eta):
    ids = np.round((x - vmin)/eta).astype(int)
    return ids[1]*numid_dim[0]+ids[0]

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
winset = xgrid[winids, :]
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
coll = PolyCollection(verts, closed=True,
                      edgecolors='k', facecolors='g', alpha=0.5)
rect_goal = patches.Rectangle((goal_real[0,0], goal_real[1,0]),
                              goal_real[0,1]-goal_real[0,0],
                              goal_real[1,1]-goal_real[1,0],
                              linewidth=1.5, edgecolor='b',fill=False)

fig, ax = plt.subplots()
ax.set_xlim(X[0,0], X[0,1])
ax.set_ylim(X[1,0], X[1,1])
ax.add_collection(coll)
ax.add_patch(rect_goal)



# # Simulation
a = 1./3.5;
B = 2.0;
H = 0.18;
W = 0.25;
lc = 8.0;
cx = 1.0/lc;
cy = 1.0/(4*lc*B*B);
aH = a+H;
H2 = H/(2.0*W*W*W);
W2 = 3*W*W;
def engineMG(t, x, u):
    dx = cx * (aH+H2*(x[0]-W)*(W2-(x[0]-W)*(x[0]-W)) - x[1]) + u[0]
    dy = cy * (x[0] - u[1]*np.sqrt(x[1]))
    return np.array([dx, dy])


rng = np.random.default_rng()
max_num_achieve_acc = 5
max_num_iteration = 500
Tsim = 50

x0 = np.array([0.5056, 0.6595])
t = 0
x = x0
q = q0
tsim = []
xsim = []
usim = []
qsim = []

while(t<Tsim and is_inside(x, X)):
    i = val_to_id(x, numid_dim, vmin, eta) #convert x to node id
    p5 = nts_ctrlr[encode3[i], :]
    p7 = ctlr[p5[2]:p5[2]+p5[0], :]
    uset = np.argwhere(p7[:, 0]==q0).squeeze()

    if(uset.size > 1):
        uid = rng.choice(uset, 1) #randomly pick one
    else:
        uid = int(uset)
    u = U[p7[uid,1], :].squeeze()

    # Integrate ode
    sol = solve_ivp(engineMG, [0, tau], x, method='RK45', args=(u,))
    tt = sol.t[-1]
    y = sol.y[:,-1]

    # Update DBA state
    q = q_prime[p5[1]*nDBA+q]

    # Save trajectories
    tsim.append(t)
    xsim.append(x)
    usim.append(u)
    d=qsim.append(q)
    # Update state
    x = y
    t += tt


xsim = np.asarray(xsim)
ax.plot(xsim[:, 0], xsim[:, 1], 'k')
ax.plot(xsim[0, 0], xsim[0, 1], marker='o', markerfacecolor='r')


plt.show()
