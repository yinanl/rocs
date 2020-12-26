from os.path import dirname, realpath
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
pypath = dirname(dirname(dirname(realpath(__file__)))) + '/python/'
sys.path.insert(1, pypath)
import utils
from odes import scara_2dbint, scara


# # Read the controller file and problem settings
dirpath = dirname(realpath(__file__))
ctlrfile = "/controller_2dbint_reach1.h5"
tau, X, U, _, _, pavings, tag, ctlr = \
    utils.read_controller_itvl_from_h5(dirpath+ctlrfile)
G = np.array([[0.498, 0.5772], [1.5739, 1.7055],
              [-0.1, 0.1], [-0.1, 0.1]])
model = scara(tau)


winset = pavings[np.argwhere(tag == 1).squeeze()]
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

avoid = pavings[np.argwhere(tag == -1).squeeze()]


# # Simulation
rng = np.random.default_rng()
Tsim = 50

x0 = np.array([0.5, -0.88, 0.05, -0.05])
t = 0
x = x0
tsim = []
xsim = []
usim = []
torqsim = []
while(t < Tsim and not utils.is_inside(x, G)):
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
    torqsim.append(torq.astype(float).T.squeeze())

    # Update state
    x = y
    t += tt

xsim = np.asarray(xsim)
usim = np.asarray(usim)
tsim = np.asarray(tsim)
torqsim = np.asarray(torqsim)

# # Save data to a .txt file
np.save(dirpath+'/traj_reach.npy', xsim)
np.save(dirpath+'/torq_reach.npy', torqsim)


# # x-y 2D plot of state space
fig, ax = plt.subplots()
ax.set_xlim(X[0, 0], X[0, 1])
ax.set_xlabel(r'$\theta_1$')
ax.set_ylim(X[1, 0], X[1, 1])
ax.set_ylabel(r'$\theta_2$')

rect_goal = patches.Rectangle((G[0, 0], G[1, 0]),
                              G[0, 1]-G[0, 0], G[1, 1]-G[1, 0],
                              linewidth=1.5, edgecolor='g', fill=False)
ax.add_patch(rect_goal)
ax.add_collection(
    utils.polycoll_interval_array(avoid[:, 0:4], True, 'k', 0.7)
)
# ax.add_collection(
#     utils.polycoll_interval_array(winset[:, 0:4], True, 'palegoldenrod', 0.7)
# )


# # Plot closed-loop trajectoryies
ax.plot(xsim[:, 0], xsim[:, 1], 'b')
ax.plot(x0[0], x0[1], marker='o', markerfacecolor='r')


# # Plot time-control curves
t_interp = np.linspace(0, tsim[-1], 200)
fu = interp1d(tsim, torqsim, axis=0, kind='previous')
torq_interp = fu(t_interp)

fig2, ax2 = plt.subplots()
ax2.plot(t_interp, torq_interp)
ax2.title.set_text('Time-Torque Curves')
ax2.legend((r'$\tau_1$', r'$\tau_2$'), loc="upper right")
ax2.set_xlabel(r'$t$(s)')
ax2.set_ylabel(r'$\tau_{1,2}$(N$\cdot$m)')

plt.show()
