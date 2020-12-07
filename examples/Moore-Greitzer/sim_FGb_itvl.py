from os.path import dirname,realpath
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
pypath = dirname(dirname(dirname(realpath(__file__)))) + '/python/'
sys.path.insert(1, pypath)
import utils
from odes import engineMG


dirpath = dirname(realpath(__file__))
w0 = "/controller_I_nba_w0.h5"
w1 = "/controller_I_nba_w1.h5"
w2 = "/controller_I_nba_w2.h5"
specfile = "/nba.txt"


# # Read the specification file
n_dba, n_props, q0, acc, q_prime = utils.read_spec_from_txt(dirpath+specfile)


# # Read the controller file and problem settings
tau, X, U, tag0, pavings0, ctlr0 = utils.read_controller_itvl_from_h5(dirpath+w0)
tau, X, U, tag1, pavings1, ctlr1 = utils.read_controller_itvl_from_h5(dirpath+w1)
tau, X, U, tag2, pavings2, ctlr2 = utils.read_controller_itvl_from_h5(dirpath+w2)
tag = [tag0, tag1, tag2]
pavings = [pavings0, pavings1, pavings2]
ctlr = [ctlr0, ctlr1, ctlr2]

winset = pavings[q0][np.argwhere(tag[q0] == 1).squeeze()]
print("\nWinning set coverage:")
wsize = np.sum((winset[:, 1]-winset[:, 0])*(winset[:, 3]-winset[:, 2]))
winper = "{:.2%}".format(wsize/((X[0, 1]-X[0, 0])*(X[1, 1]-X[1, 0])))
print(winper)

goal_real = np.array([[0.5009, 0.5069], [0.6575, 0.6635]])
obs = np.array([[0.520, 0.526], [0.658, 0.664]])


def get_propositions(x):
    if(x[0] > goal_real[0, 0] and x[0] < goal_real[0, 1] and
       x[1] > goal_real[1, 0] and x[1] < goal_real[1, 1]):
        return 1
    elif(x[0] > obs[0, 0] and x[0] < obs[0, 1] and
         x[1] > obs[1, 0] and x[1] < obs[1, 1]):
        return -1
    else:
        return 0


# # x-y 2D plot of winning set
rect_goal = patches.Rectangle((goal_real[0, 0], goal_real[1, 0]),
                              goal_real[0, 1]-goal_real[0, 0],
                              goal_real[1, 1]-goal_real[1, 0],
                              linewidth=1.5, edgecolor='g', fill=False)
rect_obs = patches.Rectangle((obs[0, 0], obs[1, 0]),
                             obs[0, 1]-obs[0, 0],
                             obs[1, 1]-obs[1, 0],
                             linewidth=1, edgecolor='k',
                             fill=True, facecolor='k')
fig, ax = plt.subplots()
ax.set_xlim(X[0, 0], X[0, 1])
ax.set_ylim(X[1, 0], X[1, 1])
ax.add_patch(rect_goal)
ax.add_patch(rect_obs)
ax.add_collection(utils.polycoll_winset_itvl(winset))


# # Simulation
rng = np.random.default_rng()
Tsim = 20

x0 = np.array([0.5343, 0.6553])
t = 0
x = x0
q = q0
tsim = []
xsim = []
usim = []
qsim = []
while(t < Tsim):
    x_id = utils.index_in_interval_array(x, pavings[q])
    if(x_id < 0):
        print("System state ")
        print(x)
        print(" is not inside the winning set.")
        break
    if(q < 0):
        print("System unsafe.")
        break

    if(any(ctlr[q][x_id, :])):
        uset = np.argwhere(ctlr[q][x_id, :]).squeeze()  # get the indices of valid input
    else:
        print("No valid control input.")
        break

    if(uset.size > 1):
        uid = rng.choice(uset, 1)  # randomly pick one
    else:
        uid = int(uset)
    u = U[uid, :].squeeze()

    # Integrate ode
    sol = solve_ivp(engineMG, [0, tau], x, method='RK45', args=(u,))
    tt = sol.t[-1]
    y = sol.y[:, -1]

    # Save trajectories
    tsim.append(t)
    xsim.append(x)
    usim.append(u)
    qsim.append(q)
    # Update state
    q = q_prime[q, get_propositions(x)]
    x = y
    t += tt


xsim = np.asarray(xsim)
ax.plot(xsim[:, 0], xsim[:, 1], 'b')
ax.plot(x0[0], x0[1], marker='o', markerfacecolor='r')


plt.show()
