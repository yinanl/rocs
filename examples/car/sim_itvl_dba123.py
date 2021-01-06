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
from odes import car
import argparse

parser = argparse.ArgumentParser(description='Indicate control scenarios.')
parser.add_argument("spec", type=str, help="indicate the case number")
args = parser.parse_args()


dirpath = dirname(realpath(__file__))
# spec = "dba2"
spec = args.spec


# # Set up workspace
theta = 3.5
nG = 3
nA = 4
G = np.zeros(shape=(3, 2, 3))
G[:, :, 0] = np.array([[1.0, 2.0], [0.5, 2.0], [-theta, theta]])
G[:, :, 1] = np.array([[0.5, 2.5], [7.5, 8.5], [-theta, theta]])
G[:, :, 2] = np.array([[7.1, 9.1], [4.6, 6.4], [-theta, theta]])
A = np.zeros(shape=(3, 2, 4))
A[:, :, 0] = np.array([[1.6, 5.7], [4.0, 5.0], [-theta, theta]])
A[:, :, 1] = np.array([[3.0, 5.0], [5.0, 8.0], [-theta, theta]])
A[:, :, 2] = np.array([[4.3, 5.7], [1.8, 4.0], [-theta, theta]])
A[:, :, 3] = np.array([[5.7, 8.5], [1.8, 2.5], [-theta, theta]])


def get_labels(x):
    if(x[0] > G[0, 0, 0] and x[0] < G[0, 1, 0] and
       x[1] > G[1, 0, 0] and x[1] < G[1, 1, 0]):
        return 4
    elif(x[0] > G[0, 0, 1] and x[0] < G[0, 1, 1] and
         x[1] > G[1, 0, 1] and x[1] < G[1, 1, 1]):
        return 2
    elif(x[0] > G[0, 0, 2] and x[0] < G[0, 1, 2] and
         x[1] > G[1, 0, 2] and x[1] < G[1, 1, 2]):
        return 1
    elif(x[0] > A[0, 0, 0] and x[0] < A[0, 1, 0] and
         x[1] > A[1, 0, 0] and x[1] < A[1, 1, 0]):
        return -1
    elif(x[0] > A[0, 0, 1] and x[0] < A[0, 1, 1] and
         x[1] > A[1, 0, 1] and x[1] < A[1, 1, 1]):
        return -1
    elif(x[0] > A[0, 0, 2] and x[0] < A[0, 1, 2] and
         x[1] > A[1, 0, 2] and x[1] < A[1, 1, 2]):
        return -1
    elif(x[0] > A[0, 0, 3] and x[0] < A[0, 1, 3] and
         x[1] > A[1, 0, 3] and x[1] < A[1, 1, 3]):
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
# x0 = np.random.rand(3)*(X[:, 1]-X[:, 0]) + X[:, 0]

xsim, usim, qsim, tsim = utils.simulate_itvl_dba_control(
    tau, Tsim, num_acc, x0, car, dba, controller, get_labels)


# # x-y 2D plot of workspace
fig, ax = plt.subplots()
ax.set_xlim(X[0, 0], X[0, 1])
ax.set_ylim(X[1, 0], X[1, 1])
FS = 16
LW = 2
for i in range(nG):
    ax.add_patch(
        patches.Rectangle((G[0, 0, i], G[1, 0, i]),
                          G[0, 1, i]-G[0, 0, i], G[1, 1, i]-G[1, 0, i],
                          facecolor='gold')
    )

for i in range(nA):
    ax.add_patch(
        patches.Rectangle((A[0, 0, i], A[1, 0, i]),
                          A[0, 1, i]-A[0, 0, i], A[1, 1, i]-A[1, 0, i],
                          facecolor='dimgray')
    )
plt.text(1.4, 1.2, '$a_1$', fontsize=FS)
plt.text(1.3, 8.0, '$a_2$', fontsize=FS)
plt.text(8.0, 5.5, '$a_3$', fontsize=FS)

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
