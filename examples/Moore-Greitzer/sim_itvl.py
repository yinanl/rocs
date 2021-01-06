import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
from os.path import dirname, realpath
pypath = dirname(dirname(dirname(realpath(__file__)))) + '/python/'
sys.path.insert(1, pypath)
import utils
from odes import MG, MG_scaled, MG_plotEP
import argparse
parser = argparse.ArgumentParser(description='Indicate control scenarios.')
parser.add_argument("case", type=int, help="indicate the case number")
args = parser.parse_args()
# from matplotlib import rc
# rc('font',**{'family':'serif','serif':['Times']})
# rc('text', usetex=True)


dirpath = dirname(realpath(__file__))
spec = "FGb"
s = 1.0


# # Set up state space
if(args.case == 1):
    # case I
    goal_real = np.array([[0.5009, 0.5069], [0.6575*s, 0.6635*s]])
    obs = np.array([[0.520, 0.526], [0.658*s, 0.664*s]])
    head = "/controller_I_"
elif(args.case == 2):
    # case II
    goal_real = np.array([[0.4489, 0.4549], [0.6483, 0.6543]])
    obs = np.array([[0.497, 0.503], [0.650, 0.656]])
    head = "/controller_II_"
else:
    raise ValueError('Wrong input value.')


def get_labels(x):
    if(x[0] > goal_real[0, 0] and x[0] < goal_real[0, 1] and
       x[1] > goal_real[1, 0] and x[1] < goal_real[1, 1]):
        return 1
    elif(x[0] > obs[0, 0] and x[0] < obs[0, 1] and
         x[1] > obs[1, 0] and x[1] < obs[1, 1]):
        return -1
    else:
        return 0


# # Read the specification file
specfile = '/' + spec + ".txt"
dba = utils.read_spec_from_txt(dirpath+specfile)


# # Read the controller files
tag = []
pavings = []
ctlr = []
for k in range(dba.n_dba):
    w = head + spec + "_w" + str(k) + ".h5"
    tau, X, U, G, A, p, t, c = utils.read_controller_itvl_from_h5(dirpath+w)
    tag.append(t)
    pavings.append(p)
    ctlr.append(c)
controller = utils.CtlrItvl(U, pavings, tag, ctlr)


# # Compute the percentage of winning set on the state space
winset = pavings[dba.q0][np.argwhere(tag[dba.q0] == 1).squeeze()]
print("\nWinning set coverage:")
wsize = np.sum((winset[:, 1]-winset[:, 0])*(winset[:, 3]-winset[:, 2]))
winper = "{:.2%}".format(wsize/((X[0, 1]-X[0, 0])*(X[1, 1]-X[1, 0])))
print(winper)


# # x-y 2D plot of winning set
fig, ax = plt.subplots()
ax.set_xlim(X[0, 0], X[0, 1])
ax.set_xlabel(r'$\Phi$', **{'fontname':'Times New Roman'})
ax.set_ylim(X[1, 0], X[1, 1])
ax.set_ylabel(r'$\Psi$', **{'fontname':'Times New Roman'})
ax.set_aspect('equal',adjustable='box')
# ax.axis('equal')
MG_plotEP(ax)
# rect_ss = patches.Rectangle((X[0, 0], X[1, 0]), X[0, 1]-X[0, 0], X[1, 1]-X[1, 0],
#                               linewidth=1.5, edgecolor='k', fill=False)
rect_goal = patches.Rectangle((goal_real[0, 0], goal_real[1, 0]),
                              goal_real[0, 1]-goal_real[0, 0],
                              goal_real[1, 1]-goal_real[1, 0],
                              linewidth=1.5, edgecolor='g', fill=False)
rect_obs = patches.Rectangle((obs[0, 0], obs[1, 0]),
                             obs[0, 1]-obs[0, 0],
                             obs[1, 1]-obs[1, 0],
                             linewidth=1, edgecolor='k',
                             fill=True, facecolor='grey')
# ax.add_patch(rect_ss)
ax.add_patch(rect_goal)
ax.add_patch(rect_obs)
ax.add_collection(
    utils.polycoll_interval_array(winset, True, 'palegoldenrod', 0.7)
)


# # Simulation
rng = np.random.default_rng()
Tsim = 20

x0 = np.array([0.5343, 0.6553*s])
t = 0
x = x0
q = dba.q0
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
    sol = solve_ivp(MG, [0, tau], x, method='RK45', args=(u,))
    tt = sol.t[-1]
    y = sol.y[:, -1]

    # Save trajectories
    tsim.append(t)
    xsim.append(x)
    usim.append(u)
    qsim.append(q)
    # Update state
    q = dba.q_prime[q, get_labels(x)]
    x = y
    t += tt


xsim = np.asarray(xsim)
ax.plot(xsim[:, 0], xsim[:, 1], 'b')
ax.plot(x0[0], x0[1], marker='o', markerfacecolor='r')


plt.show()
