import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import sys
from os.path import dirname, realpath
pypath = dirname(dirname(dirname(realpath(__file__)))) + '/python/'
sys.path.insert(1, pypath)
import utils
from odes import car


dirpath = dirname(realpath(__file__))
spec = "dba4"

# # Set up workspace
theta = 3.5
sdim = 3
udim = 2
X = np.array([[0, 10], [0, 10], [-theta, theta]])
U = np.array([[-1, 1], [-1, 1]])

nG = 6
G = np.zeros(shape=(3, 2, nG))
G[:, :, 0] = np.array([[0.0, 4.0], [6.0, 10.0], [-theta, theta]])  # c
G[:, :, 1] = np.array([[6.0, 10.0], [0.0, 4.0], [-theta, theta]])  # d
G[:, :, 2] = np.array([[1.0, 2.0], [0.5, 2.0], [-theta, theta]])  # a1
G[:, :, 3] = np.array([[0.5, 2.5], [7.5, 8.5], [-theta, theta]])  # a2
G[:, :, 4] = np.array([[7.1, 9.1], [1.9, 2.9], [-theta, theta]])  # a3
G[:, :, 5] = np.array([[3.8, 4.6], [3.1, 4.5], [-theta, theta]])  # a4
nA = 3
A = np.zeros(shape=(3, 2, nA))
A[:, :, 0] = np.array([[0.0, 3.2], [4.0, 5.0], [-theta, theta]])
A[:, :, 1] = np.array([[5.4, 6.0], [5.0, 10.0], [-theta, theta]])
A[:, :, 2] = np.array([[4.5, 5.2], [0.0, 2.5], [-theta, theta]])


def get_labels(x):
    if(x[0] > G[0, 0, 2] and x[0] < G[0, 1, 2] and
       x[1] > G[1, 0, 2] and x[1] < G[1, 1, 2]):  # a1
        return 32
    elif(x[0] > G[0, 0, 5] and x[0] < G[0, 1, 5] and
         x[1] > G[1, 0, 5] and x[1] < G[1, 1, 5]):  # a4
        return 1
    elif(x[0] > G[0, 0, 0] and x[0] < G[0, 1, 0] and
         x[1] > G[1, 0, 0] and x[1] < G[1, 1, 0]):  # c
        if(x[0] > G[0, 0, 3] and x[0] < G[0, 1, 3] and
           x[1] > G[1, 0, 3] and x[1] < G[1, 1, 3]):  # c&a2
            return 10
        else:  # c&!a2
            return 8
    elif(x[0] > G[0, 0, 1] and x[0] < G[0, 1, 1] and
         x[1] > G[1, 0, 1] and x[1] < G[1, 1, 1]):  # d
        if(x[0] > G[0, 0, 4] and x[0] < G[0, 1, 4] and
           x[1] > G[1, 0, 4] and x[1] < G[1, 1, 4]):  # d&a3
            return 20
        else:  # d&!a3
            return 4
    elif(x[0] > A[0, 0, 0] and x[0] < A[0, 1, 0] and
         x[1] > A[1, 0, 0] and x[1] < A[1, 1, 0]):
        return -1
    elif(x[0] > A[0, 0, 1] and x[0] < A[0, 1, 1] and
         x[1] > A[1, 0, 1] and x[1] < A[1, 1, 1]):
        return -1
    elif(x[0] > A[0, 0, 2] and x[0] < A[0, 1, 2] and
         x[1] > A[1, 0, 2] and x[1] < A[1, 1, 2]):
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
    w = "/controller_" + spec + "_w" + str(k) + ".h5"
    tau, X, U, _, _, p, t, c = utils.read_controller_itvl_from_h5(dirpath+w)
    tag.append(t)
    pavings.append(p)
    ctlr.append(c)
controller = utils.CtlrItvl(U, pavings, tag, ctlr)


# # Compute the percentage of winning set on the state space
winset = pavings[dba.q0][np.argwhere(tag[dba.q0] == 1).squeeze()]
print("\nWinning set coverage:")
wsize = np.sum((winset[:, 1]-winset[:, 0])
               * (winset[:, 3]-winset[:, 2])
               * (winset[:, 5]-winset[:, 4]))
winper = "{:.2%}".format(
    wsize/((X[0, 1]-X[0, 0])*(X[1, 1]-X[1, 0])*(X[2, 1]-X[2, 0]))
)
print(winper)


# # Simulation
Tsim = 50
num_acc = 3
x0 = np.array([3.0, 2.0, np.pi/2.])
# x0 = np.array([3.7, 5.3, -np.pi/2.])
xsim, usim, qsim, tsim = utils.simulate_itvl_dba_control(
    tau, Tsim, num_acc, x0, car, dba, controller, get_labels)


# # x-y 2D plot of winning set
FS = 16
LW = 2

fig, ax = plt.subplots()
ax.set_xlim(X[0, 0], X[0, 1])
ax.set_ylim(X[1, 0], X[1, 1])

for i in range(nG):
    if(i < 2):
        clr = 'aquamarine'
    else:
        clr = 'gold'
    ax.add_patch(
        patches.Rectangle((G[0, 0, i], G[1, 0, i]),
                          G[0, 1, i]-G[0, 0, i], G[1, 1, i]-G[1, 0, i],
                          facecolor=clr)
    )

for i in range(nA):
    ax.add_patch(
        patches.Rectangle((A[0, 0, i], A[1, 0, i]),
                          A[0, 1, i]-A[0, 0, i], A[1, 1, i]-A[1, 0, i],
                          facecolor='dimgray')
    )

plt.text(1.4, 1.2, '$a_1$', fontsize=FS)
plt.text(1.3, 8.0, '$a_2$', fontsize=FS)
plt.text((G[0, 0, 4]+G[0, 1, 4])/2,(G[1, 1, 4]+G[1, 0, 4])/2, '$a_3$',fontsize=FS)
plt.text((G[0, 0, 5]+G[0, 1, 5])/2,(G[1, 1, 5]+G[1, 0, 5])/2, '$a_4$',fontsize=FS)
plt.text(3.5, 6.5, '$c$', fontsize=FS)
plt.text(6.5, 3.5, '$d$', fontsize=FS)

# ax.add_collection(
#     utils.polycoll_interval_array(winset, True, 'palegoldenrod', 0.7)
# )


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
fu = interp1d(tsim, usim, axis=0, kind='previous')
u_interp = fu(t_interp)
fq = interp1d(tsim, qsim, axis=0, kind='previous')
q_interp = fq(t_interp)

axs[0].plot(t_interp, u_interp)
axs[0].title.set_text('Time-Control Curves')
axs[0].legend(('$v$', '$\omega$'), loc="upper right")
axs[1].plot(t_interp, q_interp)
axs[1].title.set_text('Time-Automaton State Curve')
plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                    wspace=None, hspace=0.5)


plt.show()
