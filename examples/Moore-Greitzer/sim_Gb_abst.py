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
filename = "/controller_Gb_abst_0.00018.h5"
specfile = "/Gb.txt"


# # Read the specification file
n_dba, _, q0, acc, _ = utils.read_spec_from_txt(dirpath+specfile)


# # Read controller file and problem settings
tau, X, U, eta, xgrid, goalset, winids, ctlr, encode3, nts_ctrlr, q_prime = \
    utils.read_controller_abst_from_h5(dirpath+filename)
winset = xgrid[winids, :]
print("\nWinning set coverage:")
winper = "{:.2%}".format(winids.size/xgrid.shape[0])
print(winper)

goal_real = np.array([[0.5009, 0.5069], [0.6575, 0.6635]])


# # x-y 2D plot of winning set
rect_goal = patches.Rectangle((goal_real[0, 0], goal_real[1, 0]),
                              goal_real[0, 1]-goal_real[0, 0],
                              goal_real[1, 1]-goal_real[1, 0],
                              linewidth=1.5, edgecolor='b', fill=False)
fig, ax = plt.subplots()
ax.set_xlim(X[0, 0], X[0, 1])
ax.set_ylim(X[1, 0], X[1, 1])
ax.add_collection(utils.polycoll_winset_abst(winset, eta))
ax.add_patch(rect_goal)


# # Simulation
rng = np.random.default_rng()
Tsim = 50

x0 = np.array([0.5056, 0.6595])
t = 0
x = x0
q = q0
tsim = []
xsim = []
usim = []
qsim = []
while(t < Tsim):
    i = utils.index_in_grid(x, xgrid)  # convert x to node id

    p5 = nts_ctrlr[encode3[i], :]
    p7 = ctlr[p5[2]:p5[2]+p5[0], :]
    uset = np.argwhere(p7[:, 0] == q0).squeeze()

    if(uset.size > 1):
        uid = rng.choice(uset, 1)  # randomly pick one
    else:
        uid = int(uset)
    u = U[p7[uid, 1], :].squeeze()

    # Integrate ode
    sol = solve_ivp(MG, [0, tau], x, method='RK45', args=(u,))
    tt = sol.t[-1]
    y = sol.y[:, -1]

    # Update DBA state
    q = q_prime[p5[1]*n_dba+q]

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
