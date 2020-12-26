from os.path import dirname, realpath
import numpy as np
from functools import reduce
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
pypath = dirname(dirname(dirname(realpath(__file__)))) + '/python/'
sys.path.insert(1, pypath)
import utils
from odes import dcdc


# # Read the controller file and problem settings
dirpath = dirname(realpath(__file__))
# ctlrfile = "/controller_dcdcInv.h5"
ctlrfile = "/controller_dcdcCoBuchi.h5"
tau, X, U, G, _, pavings, tag, ctlr = \
    utils.read_controller_itvl_from_h5(dirpath+ctlrfile)
G = G.squeeze()

winset = pavings[np.argwhere(tag == 1).squeeze()]
print("\nWinning set coverage:")
wsize = np.sum((winset[:, 1]-winset[:, 0])*(winset[:, 3]-winset[:, 2]))
winper = "{:.2%}".format(
    wsize/((X[0, 1]-X[0, 0])*(X[1, 1]-X[1, 0]))
)
print(winper)

avoid = pavings[np.argwhere(tag == -1).squeeze()]


# # Simulation
Tsim = 50
# x0 = np.array([1.2, 1.12])
x0 = np.array([0.7, 5.4/5])

t = 0
x = x0
tsim = []
xsim = []
usim = []
rng = np.random.default_rng()
while(t < Tsim):
    x_id = utils.index_in_interval_array(x, pavings)
    if(x_id < 0):
        print("System state ")
        print(x)
        print(" is not inside the winning set.")
        break

    if(any(ctlr[x_id, :])):
        uset = np.argwhere(ctlr[x_id, :]).squeeze()  # get the indices of valid input
    else:
        print("No valid control input.")
        break

    if(uset.size > 1):
        uid = rng.choice(uset, 1)  # randomly pick one
    else:
        uid = int(uset)
    u = U[uid, :].squeeze()

    # Integrate ode
    y = dcdc(tau, np.atleast_2d(x).T, u)  # another transpose 1D: x[..., None]

    # Save trajectories
    tsim.append(t)
    xsim.append(x)
    usim.append(u)
    # Update state
    x = y.T.squeeze()
    t += tau

xsim = np.asarray(xsim)
usim = np.asarray(usim)
tsim = np.asarray(tsim)


# # x-y 2D plot of state space
fig, ax = plt.subplots()
ax.set_xlim(X[0, 0], X[0, 1])
ax.set_ylim(X[1, 0], X[1, 1])

rect_goal = patches.Rectangle((G[0, 0], G[1, 0]), G[0, 1]-G[0, 0], G[1, 1]-G[1, 0],
                              linewidth=1.5, edgecolor='g', fill=False)
ax.add_patch(rect_goal)
ax.add_collection(
    utils.polycoll_interval_array(avoid[:, 0:4], True, 'k', 0.7)
)
ax.add_collection(
    utils.polycoll_interval_array(winset[:, 0:4], True, 'palegoldenrod', 0.7)
)


# # Plot closed-loop trajectoryies
ax.plot(xsim[:, 0], xsim[:, 1], 'b')
ax.plot(x0[0], x0[1], marker='o', markerfacecolor='r')


plt.show()
