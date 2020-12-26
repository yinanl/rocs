import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import sys
from os.path import dirname, realpath
pypath = dirname(dirname(dirname(realpath(__file__)))) + '/python/'
sys.path.insert(1, pypath)
import utils
from odes import scara_2dbint, scara


# # Read the specification file
dirpath = dirname(realpath(__file__))
spec = "gb2"
specfile = '/' + spec + ".txt"
dba = utils.read_spec_from_txt(dirpath+specfile)


# # Read the controller files
tag = []
pavings = []
ctlr = []
for k in range(dba.n_dba):
    w = "/controller_" + spec + "_w" + str(k) + ".h5"
    tau, X, U, _, _, p, t, c = utils.read_controller_itvl_from_h5(dirpath+w)
    tag.append(t)
    pavings.append(p)
    ctlr.append(c)
controller = utils.CtlrItvl(U, pavings, tag, ctlr)
model = scara(tau)


# # Compute the percentage of winning set on the state space
winset = pavings[dba.q0][np.argwhere(tag[dba.q0] == 1).squeeze()]
print("\nWinning set coverage:")
wsize = np.sum((winset[:, 1]-winset[:, 0])
               * (winset[:, 3]-winset[:, 2])
               * (winset[:, 5]-winset[:, 4])
               * (winset[:, 7]-winset[:, 6]))
winper = "{:.2%}".format(
    wsize/((X[0, 1]-X[0, 0])*(X[1, 1]-X[1, 0])
           * (X[2, 1]-X[2, 0])*(X[3, 1]-X[3, 0]))
)
print(winper)


# # Set up workspace
nG = 2
G = np.zeros(shape=(4, 2, 2))
G[:, :, 0] = np.array([[0.4980, 0.5772], [1.5739, 1.7055],
                       [-0.1, 0.1], [-0.1, 0.1]])
G[:, :, 1] = np.array([[0.4903, 0.6069], [-0.9363, -0.8363],
                       [-0.1, 0.1], [-0.1, 0.1]])

A = pavings[dba.q0][np.argwhere(tag[dba.q0] == -1).squeeze(), :]


def get_labels(x, G, A):
    if(x[0] > G[0, 0, 0] and x[0] < G[0, 1, 0] and
       x[1] > G[1, 0, 0] and x[1] < G[1, 1, 0]):
        return 1
    elif(x[0] > G[0, 0, 1] and x[0] < G[0, 1, 1] and
         x[1] > G[1, 0, 1] and x[1] < G[1, 1, 1]):
        return 2
    elif(utils.index_in_interval_array(x, A).size > 0):
        return -1
    else:
        return 0


# # Simulation
rng = np.random.default_rng()
Tsim = 50

x0 = np.array([0.0, 0.0, 0.0, 0.0])
t = 0
x = x0
q = dba.q0
tsim = []
xsim = []
usim = []
qsim = []
torqsim = []
lsim = []
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

    # Calculate torque
    torq = model.compute_torque(np.atleast_2d(u).T, np.atleast_2d(x[2:4]).T,
                                np.atleast_2d(x[0:2]).T)

    # Integrate ode
    sol = solve_ivp(scara_2dbint, [0, tau], x, method='RK45', args=(u,))
    tt = sol.t[-1]
    y = sol.y[:, -1]

    # Save trajectories
    tsim.append(t)
    xsim.append(x)
    usim.append(u)
    qsim.append(q)
    torqsim.append(torq.astype(float).T.squeeze())
    lsim.append(get_labels(x, G, A))
    # Update state
    q = dba.q_prime[q, get_labels(x, G, A)]
    x = y
    t += tt

xsim = np.asarray(xsim)
usim = np.asarray(usim)
qsim = np.asarray(qsim)
tsim = np.asarray(tsim)
torqsim = np.asarray(torqsim)
lsim = np.asarray(lsim)


# # Save data to a .txt file
np.save(dirpath+'/traj_gb2.npy', xsim)
np.save(dirpath+'/torq_gb2.npy', torqsim)


# # x-y 2D plot of winning set
fig, ax = plt.subplots()
ax.set_xlim(X[0, 0], X[0, 1])
ax.set_xlabel(r'$\theta_1$')
ax.set_ylim(X[1, 0], X[1, 1])
ax.set_ylabel(r'$\theta_2$')

for i in range(nG):
    ax.add_patch(
        patches.Rectangle((G[0, 0, i], G[1, 0, i]),
                          G[0, 1, i]-G[0, 0, i], G[1, 1, i]-G[1, 0, i],
                          linewidth=1.5, edgecolor='g', fill=False)
    )
ax.add_collection(
    utils.polycoll_interval_array(A[:, 0:4], True, 'k', 0.7)
)
ax.add_collection(
    utils.polycoll_interval_array(winset[:, 0:4], True, 'palegoldenrod', 0.7)
)


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


figs, axs = plt.subplots(2)
t_interp = np.linspace(0, tsim[-1], 200)
fu = interp1d(tsim, torqsim, axis=0, kind='previous')
torq_interp = fu(t_interp)
fq = interp1d(tsim, qsim, axis=0, kind='previous')
q_interp = fq(t_interp)

axs[0].plot(t_interp, torq_interp)
axs[0].title.set_text('Time-Control Curves')
axs[0].legend((r'$\tau_1$', r'$\tau_2$'), loc="upper right")
axs[0].set_xlabel(r'$t$(s)')
axs[0].set_ylabel(r'$\tau_{1,2}$(N$\cdot$m)')
axs[1].plot(t_interp, q_interp)
axs[1].title.set_text('Time-Automaton State Curve')
axs[1].set_xlabel(r'$t$(s)')
axs[1].set_ylabel(r'$q$')
plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                    wspace=None, hspace=0.5)

plt.show()
