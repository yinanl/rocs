from os.path import dirname,realpath
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
pypath = dirname(dirname(dirname(realpath(__file__)))) + '/python/'
sys.path.insert(1, pypath)
import utils
from odes import MG


dirpath = dirname(realpath(__file__))
ctlrfile = "/controller_abstI_0.00009.h5"
specfile = "/nba.txt"


# # Set up state space
goal_real = np.array([[0.5009, 0.5069], [0.6575, 0.6635]])
obs = np.array([[0.520, 0.526], [0.658, 0.664]])


# # Read the specification file
dba = utils.read_spec_from_txt(dirpath+specfile)


# # Read controller file and problem settings
tau, X, eta, _, winids, controller = \
    utils.read_controller_abst_from_h5(dirpath+ctlrfile)


# # Compute the percentage of winning set on the state space
winset = controller.xgrid[winids, :]
print("\nWinning set coverage:")
winper = "{:.2%}".format(winids.size/controller.xgrid.shape[0])
print(winper)



# # x-y 2D plot of state space
fig, ax = plt.subplots()
ax.set_xlim(X[0, 0], X[0, 1])
ax.set_ylim(X[1, 0], X[1, 1])
rect_goal = patches.Rectangle((goal_real[0, 0], goal_real[1, 0]),
                              goal_real[0, 1]-goal_real[0, 0],
                              goal_real[1, 1]-goal_real[1, 0],
                              linewidth=1.5, edgecolor='g', fill=False)
rect_obs = patches.Rectangle((obs[0, 0], obs[1, 0]),
                             obs[0, 1]-obs[0, 0],
                             obs[1, 1]-obs[1, 0],
                             linewidth=1, edgecolor='k',
                             fill=True, facecolor='k')
ax.add_patch(rect_goal)
ax.add_patch(rect_obs)
# ax.add_collection(
#     utils.polycoll_grid_array(winset, eta, True, 'palegoldenrod', 0.7)
# )


# # Simulation
Tsim = 20
x0 = np.array([0.5343, 0.6553])

rng = np.random.default_rng()
t = 0
x = x0
q = dba.q0
tsim = []
xsim = []
usim = []
qsim = []

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
    sol = solve_ivp(MG, [0, tau], x, method='RK45', args=(u,))
    tt = sol.t[-1]
    y = sol.y[:, -1]

    # Update DBA state
    q = controller.q_prime[p5[1]*dba.n_dba+q]  # p5[1] is the label/proposition of current x

    # Save trajectories
    tsim.append(t)
    xsim.append(x)
    usim.append(u)
    qsim.append(q)
    # Update state
    x = y
    t += tt


xsim = np.asarray(xsim)
ax.plot(xsim[:, 0], xsim[:, 1], 'b')
ax.plot(xsim[0, 0], xsim[0, 1], marker='o', markerfacecolor='r')


plt.show()
