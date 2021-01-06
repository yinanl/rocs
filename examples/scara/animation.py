import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines
import matplotlib.animation as animation
import sys
from os.path import dirname, realpath
pypath = dirname(dirname(dirname(realpath(__file__)))) + '/python/'
sys.path.insert(1, pypath)
from odes import scara


# # Load simulated trajectories
dirpath = dirname(realpath(__file__))
statefile = "traj_gb2"
torqfile = "torq_gb2"
sfile = dirpath + '/' + statefile + ".npy"
tfile = dirpath + '/' + torqfile + ".npy"

thetas = np.load(sfile)
torqs = np.load(tfile)

# # Joint space to operation space
model = scara(0.1)
xy2 = model.theta2xy(thetas[:, 0], thetas[:, 1])
x1 = model.l1 * np.cos(thetas[:, 0])
y1 = model.l2 * np.sin(thetas[:, 0])
x2 = xy2[0, :]
y2 = xy2[1, :]


# # Setup workspace
FS = 12

fig, ax = plt.subplots()
ax.set_xlim(0, 1.2*(model.l1+model.l2))
ax.set_xlabel(r'$x$')
ax.set_ylim(0, model.l1+model.l2)
ax.set_ylabel(r'$y$')
ax.title

# The bar obstacle
H = 0.8*model.l1
r = 0.5*model.l1
# bar = plt.plot([0, r], [H, H], linewidth=10, color='k')
bar = patches.Rectangle((0, H), r, 0.01, facecolor='tab:gray',
                        hatch='/', zorder=0)
ax.add_patch(bar)

# Goals
nG = 2
G = np.zeros(shape=(2, 2, nG))
G[:, :, 0] = np.array([[0.0277, 0.0597], [0.1852, 0.2134]])
G[:, :, 1] = np.array([[0.2585, 0.2784], [0.0059, 0.0514]])
for i in range(nG):
    ax.add_patch(
        patches.Rectangle((G[0, 0, i], G[1, 0, i]),
                          G[0, 1, i]-G[0, 0, i], G[1, 1, i]-G[1, 0, i],
                          linewidth=1.5, facecolor='yellowgreen', fill=True,
                          zorder=0)
    )
plt.text((G[0, 0, 0]+G[0, 1, 0])*0.4, (G[1, 1, 0]+G[1, 0, 0])/2,
         r'$g_1$', fontsize=FS)
plt.text((G[0, 0, 1]+G[0, 1, 1])*0.49, (G[1, 1, 1]+G[1, 0, 1])/2,
         r'$g_2$', fontsize=FS)

# arm1 = lines.Line2D([0, x1[0]], [0, y1[0]],
#                     linewidth=3, color='k', alpha=0, zorder=1)
# arm2 = lines.Line2D([x1[0], x2[0]], [y1[0], y2[0]],
#                     linewidth=3, color='k', alpha=0, zorder=1)
# joint1 = patches.Circle((0, 0), radius=0.005, color='k',
#                         fill=True, alpha=1, zorder=2)
# joint2 = patches.Circle((x1[0], y1[0]), radius=0.005, color='k',
#                         fill=True, alpha=0, zorder=2)
# end = patches.Circle((x2[0], y2[0]), radius=0.005, color='tab:orange',
#                      fill=True, alpha=0, zorder=2)
i = 93
arm1 = lines.Line2D([0, x1[i]], [0, y1[i]],
                    linewidth=3, color='k', alpha=1, zorder=1)
arm2 = lines.Line2D([x1[i], x2[i]], [y1[i], y2[i]],
                    linewidth=3, color='k', alpha=1, zorder=1)
joint1 = patches.Circle((0, 0), radius=0.005, color='k',
                        fill=True, alpha=1, zorder=2)
joint2 = patches.Circle((x1[i], y1[i]), radius=0.005, color='k',
                        fill=True, alpha=1, zorder=2)
end = patches.Circle((x2[i], y2[i]), radius=0.005, color='tab:orange',
                     fill=True, alpha=1, zorder=2)
ax.add_patch(joint1)
ax.add_patch(joint2)
ax.add_patch(end)
ax.add_artist(arm1)
ax.add_artist(arm2)


# # Animation
torque_text = ax.text(0.05, 0.95, '', transform=ax.transAxes)
torque_template = 'torques=%.3f,%.3f'


def animate(i):
    arm1.set_data([0, x1[i]], [0, y1[i]])
    arm2.set_data([x1[i], x2[i]], [y1[i], y2[i]])
    joint2.center = (x1[i], y1[i])
    end.center = (x2[i], y2[i])
    arm1.set_alpha(1)
    arm2.set_alpha(1)
    joint2.set_alpha(1)
    end.set_alpha(1)
    joint2.set_zorder(10)
    end.set_zorder(10)
    torque_text.set_text(torque_template % (torqs[i, 0], torqs[i, 1]))
    return joint2, end, torque_text, arm1, arm2

ani = animation.FuncAnimation(fig, animate, x1.size,
                              interval=0.1*500, blit=True)


# # End-effector trajectory
ax.plot(x2, y2, color='peru')

plt.savefig(dirpath+'/fig_traj-gb2-os.png')

plt.show()
