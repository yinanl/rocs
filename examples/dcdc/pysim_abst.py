import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import sys
from os.path import dirname, realpath
pypath = dirname(dirname(dirname(realpath(__file__)))) + '/python/'
sys.path.insert(1, pypath)
import utils
from odes import dcdc
import argparse

parser = argparse.ArgumentParser(description='Indicate control scenarios.')
parser.add_argument("spec", type=str, help="indicate the case number")
args = parser.parse_args()
spec = args.spec

dirpath = dirname(realpath(__file__))


G = np.array([[1.15, 1.55], [1.09, 1.17]])

# # Read the controller files
if(spec == 'inv'):
        ctlrfile = '/controller_abst_Gb.h5'
        specfile = '/Gb.txt'
elif(spec == 'rs'):
        ctlrfile = '/controller_abst_FGb.h5'
        specfile = '/FGb.txt'
else:
    raise ValueError('Wrong input value.')

# # Read the specification file
dba = utils.read_spec_from_txt(dirpath+specfile)

tau, X, eta, _, winids, controller = \
    utils.read_controller_abst_from_h5(dirpath+ctlrfile)
winset = controller.xgrid[winids, :]
print("\nWinning set coverage:")
winper = "{:.2%}".format(winids.size/controller.xgrid.shape[0])
print(winper)


# # Simulation
Tsim = 50
if(spec == 'inv'):
    x0 = np.array([1.2, 1.12])
elif(spec == 'rs'):
    x0 = np.array([0.7, 5.4/5])
else:
    raise ValueError('Wrong input value.')


t = 0
x = x0
q = dba.q0
tsim = []
xsim = []
usim = []
qsim = []
rng = np.random.default_rng()
while(t < Tsim):
    i = utils.index_in_grid(x, controller.xgrid)  # convert x to node id

    p5 = controller.nts_ctrlr[controller.encode3[i], :]
    p7 = controller.ctlr[p5[2]:p5[2]+p5[0], :]
    uset = np.argwhere(p7[:, 0] == q).squeeze()

    if(uset.size > 1):
        uid = rng.choice(uset, 1)  # randomly pick one
    else:
        uid = int(uset)
    u = controller.ugrid[p7[uid, 1], :].squeeze()

    # Integrate ode
    y = dcdc(tau, np.atleast_2d(x).T, u)  # another transpose 1D: x[..., None]

    # Update DBA state
    q = controller.q_prime[p5[1]*dba.n_dba+q]  # p5[1] is the label/proposition of current x

    # Save trajectories
    tsim.append(t)
    xsim.append(x)
    usim.append(u)
    qsim.append(q)
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
    utils.polycoll_grid_array(winset, eta, True, 'palegoldenrod', 0.7)
)


# # Plot closed-loop trajectoryies
ax.plot(xsim[:, 0], xsim[:, 1], 'b')
ax.plot(x0[0], x0[1], marker='o', markerfacecolor='r')


plt.show()
