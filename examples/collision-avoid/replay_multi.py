import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
import re
import numpy as np
import math
import h5py
import argparse
import sys
from os.path import dirname, realpath
pypath = dirname(dirname(dirname(realpath(__file__)))) + '/python/'
sys.path.insert(1, pypath)
import utils

parser = argparse.ArgumentParser(description='Indicate control scenarios.')
# parser.add_argument('case', type=str, help='Indicate the case number')
parser.add_argument('tsim', type=float, help='Indicate the simulation time step')
args = parser.parse_args()


dirpath = dirname(realpath(__file__))


# # Display workspace
goalclr = 'lightgreen'
gclr = 'olivedrab'
lfont = {'fontname':'Consolas'}
tfont = {'fontname':'Times New Roman'}

fig = plt.figure()
ax = plt.axes()
ax.set_xlim(0,10)
ax.set_ylim(0,10)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_title('2D Collision-free Trajectory', **tfont)
# ax.set_xticks(np.arange(0, 10, 0.2)) # Show grids
# ax.set_yticks(np.arange(0, 10, 0.2))
# plt.grid()


transfile = '/abstfull_0.2-0.2-0.2.h5'
labelfile = '/labels_dba1_abstfull_0.2-0.2-0.2.h5'
xgrid = np.array([])
eta = np.array([])
labels = np.array([])
obs = np.array([])
with h5py.File(dirpath+transfile, 'r') as ft,\
     h5py.File(dirpath+labelfile, 'r') as fl:
    eta = ft['eta'][...]
    xgrid = ft['xgrid'][...]
    obs = ft['obs'][...]
    labels = fl['labels'][...]

# oset = xgrid[obs, :]
# ax.add_collection(
#     utils.polycoll_grid_array(oset, eta, True, 'gray', 0.7)
# )

obstacles = [[0.0, 0.5, 5.0, 6.0],
             [2.4, 2.6, 0.0, 3.2],
             [3.9, 4.1, 9.0, 10.0], #[3.9, 4.1, 8.0, 10.0]
             [5.9, 6.1, 0.0, 0.6],
             [5.9, 6.1, 3.8, 5.1], # [5.9, 6.1, 3.8, 6.1],
             [6.1, 10.0, 4.9, 5.1]] # [6.1, 10.0, 5.9, 6.1]
rects_obs = [patches.Rectangle((obstacles[i][0], obstacles[i][2]),
                               obstacles[i][1]-obstacles[i][0],
                               obstacles[i][3]-obstacles[i][2],
                               linewidth=1,edgecolor='k',facecolor='k')
              for i in range(len(obstacles))]

goals = [[0.5, 2.0, 7.5, 9.5],  # a
         [7.5, 9.5, 0.8, 3.0]]  # d
rects_gs = [patches.Rectangle((goals[i][0], goals[i][2]),
                              goals[i][1]-goals[i][0],
                              goals[i][3]-goals[i][2],
                              linewidth=1,edgecolor=goalclr,#fill=False
                              facecolor=goalclr, alpha=0.7)
            for i in range(len(goals))]
circ_gs = patches.Circle([8.0, 8.0], radius=0.8, linewidth=1,
                         edgecolor=goalclr, #fill=False
                         facecolor=goalclr, alpha=0.7)

# gset = xgrid[np.where(labels>0),:].squeeze()
# ax.add_collection(
#     utils.polycoll_grid_array(gset, eta, True, gclr, 1)
# )


for rect in rects_obs:
    ax.add_patch(rect)

for rect in rects_gs:
    ax.add_patch(rect)
ax.add_patch(circ_gs)

ax.text((goals[0][0]+goals[0][1])/2.0-0.5, (goals[0][2]+goals[0][3])/2.0,
        "amount", **lfont)
ax.text((goals[1][0]+goals[1][1])/2.0-0.4, (goals[1][2]+goals[1][3])/2.0,
        "drop", **lfont)
ax.text(7.6, 7.8, "count", **lfont)


# # load simulated trajectories
egofile = dirpath+'/traj_closedloop_ca_4.txt'
otherfile1 = dirpath + '/traj_other_robot_4_1.txt'
otherfile2 = dirpath + '/traj_other_robot_4_2.txt'
otherfile3 = dirpath + '/traj_other_robot_4_3.txt'
x = []
u = []
xo1 = []
xo2 = []
xo3 = []
with open(egofile, "r") as ego:
    lines = ego.readlines()
    for line in lines:
        line = line.strip()
        xu = re.split(';', line)
        x.append([float(e) for e in re.split(',', xu[0])])
        u.append([float(e) for e in re.split(',|\n', xu[1])])
with open(otherfile1, "r") as other1, \
     open(otherfile2, "r") as other2, \
     open(otherfile3, "r") as other3:
    olines = other1.readlines()
    for line in olines:
        line = line.strip()
        xu = re.split(';', line)
        xo1.append([float(e) for e in re.split(',', xu[0])])
        # uo.append([float(e) for e in re.split(',|\n', xu[1])])
    olines = other2.readlines()
    for line in olines:
        line = line.strip()
        xu = re.split(';', line)
        xo2.append([float(e) for e in re.split(',', xu[0])])
        # uo.append([float(e) for e in re.split(',|\n', xu[1])])
    olines = other3.readlines()
    for line in olines:
        line = line.strip()
        xu = re.split(';', line)
        xo3.append([float(e) for e in re.split(',', xu[0])])
        # uo.append([float(e) for e in re.split(',|\n', xu[1])])


x = np.array(x)
u = np.array(u)
xo1 = np.array(xo1)
xo2 = np.array(xo2)
xo3 = np.array(xo3)
# print(x)
# print(u)
# print(xo)


# # animation
tau = args.tsim
width = 0.6
height = 0.5  # RobotCar size
dmin = 0.8 # minimum distance to the other robot

rdiag = np.sqrt(width**2+height**2)/2.0
a0 = np.arctan2(height, width)+np.pi
lowleft_ego = (x[0, 0]+rdiag*np.cos(a0+x[0, 2]),
               x[0, 1]+rdiag*np.sin(a0+x[0, 2]))
rec_ego = patches.Rectangle(lowleft_ego, width, height,
                            angle=x[0, 2]*180/math.pi,
                            fill = True, ec='b', fc='b', alpha=0)
r_ego = patches.Circle((x[0,0], x[0,1]), radius=dmin/2.0,
                       color='r', fill=False, alpha=0)
ax.add_patch(rec_ego)
ax.add_patch(r_ego)

lowleft_other1 = (xo1[0, 0]+rdiag*np.cos(a0+xo1[0, 2]),
                 xo1[0, 1]+rdiag*np.sin(a0+xo1[0, 2]))
rec_other1 = patches.Rectangle(lowleft_other1, width, height,
                              angle=xo1[0, 2]*180/np.pi,
                              fill = True, ec='grey', fc='grey', alpha=0)
ax.add_patch(rec_other1)
lowleft_other2 = (xo2[0, 0]+rdiag*np.cos(a0+xo2[0, 2]),
                 xo2[0, 1]+rdiag*np.sin(a0+xo2[0, 2]))
rec_other2 = patches.Rectangle(lowleft_other2, width, height,
                              angle=xo2[0, 2]*180/np.pi,
                              fill = True, ec='grey', fc='grey', alpha=0)
ax.add_patch(rec_other2)
lowleft_other3 = (xo3[0, 0]+rdiag*np.cos(a0+xo3[0, 2]),
                 xo3[0, 1]+rdiag*np.sin(a0+xo3[0, 2]))
rec_other3 = patches.Rectangle(lowleft_other3, width, height,
                              angle=xo3[0, 2]*180/np.pi,
                              fill = True, ec='grey', fc='grey', alpha=0)
ax.add_patch(rec_other3)

time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes)


def animate(i):
    time_text.set_text(time_template % (i * tau))
    rec_ego.xy = (x[i, 0]+rdiag*np.cos(a0+x[i, 2]),
                  x[i, 1]+rdiag*np.sin(a0+x[i, 2]))
    rec_ego.angle = x[i, 2]*180/np.pi
    r_ego.center = (x[i,0], x[i,1])
    rec_ego.set_alpha(0.7)
    r_ego.set_alpha(0.7)

    rec_other1.xy = (xo1[i, 0]+rdiag*np.cos(a0+xo1[i, 2]),
                    xo1[i, 1]+rdiag*np.sin(a0+xo1[i, 2]))
    rec_other1.angle = xo1[i, 2]*180/np.pi
    rec_other1.set_alpha(0.7)
    rec_other2.xy = (xo2[i, 0]+rdiag*np.cos(a0+xo2[i, 2]),
                    xo2[i, 1]+rdiag*np.sin(a0+xo2[i, 2]))
    rec_other2.angle = xo2[i, 2]*180/np.pi
    rec_other2.set_alpha(0.7)
    rec_other3.xy = (xo3[i, 0]+rdiag*np.cos(a0+xo3[i, 2]),
                    xo3[i, 1]+rdiag*np.sin(a0+xo3[i, 2]))
    rec_other3.angle = xo3[i, 2]*180/np.pi
    rec_other3.set_alpha(0.7)
    return rec_ego, rec_other1, rec_other2, rec_other3, r_ego, time_text

ani = animation.FuncAnimation(fig, animate, range(len(x)),
                              interval=tau*500, blit=True)

# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save(dirpath+'case4.mp4', writer=writer)


# # Plot nominal trajectory
egonom = dirpath + "/traj_closedloop_nominal.txt"
xnom = []
unom = []
with open(egonom, "r") as ego:
    lines = ego.readlines()
    for line in lines:
        line = line.strip()
        xu = re.split(';', line)
        xnom.append([float(e) for e in re.split(',', xu[0])])
        unom.append([float(e) for e in re.split(',|\n', xu[1])])
xnom = np.array(xnom)
unom = np.array(unom)
ax.plot(xnom[:, 0], xnom[:, 1], 'tab:orange')

# # Plot CA trajectory
ax.plot(x[:, 0], x[:, 1])
ax.plot(x[0,0], x[0,1], marker='^', markerfacecolor='r')
ax.plot(x[-1,0], x[-1,1], marker='v', markerfacecolor='g')
# ax.plot(x[76,0], x[76,1], marker='o', markerfacecolor='y')

ax.plot(xo1[xo1[:,0]>0, 0], xo1[xo1[:,0]>0, 1], 'k')
ax.plot(xo1[-1,0], xo1[-1,1], marker='v', markerfacecolor='g')
ax.plot(xo2[xo2[:,0]>0, 0], xo2[xo2[:,0]>0, 1], 'k')
ax.plot(xo2[-1,0], xo2[-1,1], marker='v', markerfacecolor='g')
ax.plot(xo3[xo3[:,0]>0, 0], xo3[xo3[:,0]>0, 1], 'k')
ax.plot(xo3[-1,0], xo3[-1,1], marker='v', markerfacecolor='g')

# # # Plot robots
# j = 150
# ll_ego = (x[j, 0]+rdiag*np.cos(a0+x[j, 2]),
#                x[j, 1]+rdiag*np.sin(a0+x[j, 2]))
# show_ego = patches.Rectangle(ll_ego, width, height,
#                             angle=x[j, 2]*180/math.pi,
#                             fill = True, ec='b', fc='b', alpha=0.7)
# rshow_ego = patches.Circle((x[j,0], x[j,1]), radius=dmin/2.0,
#                        color='r', fill=False, alpha=0.7)

# ll_other1 = (xo1[j, 0]+rdiag*np.cos(a0+xo1[j, 2]),
#                  xo1[j, 1]+rdiag*np.sin(a0+xo1[j, 2]))
# show_other1 = patches.Rectangle(ll_other1, width, height,
#                               angle=xo1[j, 2]*180/math.pi,
#                               fill = True, ec='grey', fc='grey', alpha=0.7)
# ll_other2 = (xo2[j, 0]+rdiag*np.cos(a0+xo2[j, 2]),
#                  xo2[j, 1]+rdiag*np.sin(a0+xo2[j, 2]))
# show_other2 = patches.Rectangle(ll_other2, width, height,
#                               angle=xo2[j, 2]*180/math.pi,
#                               fill = True, ec='grey', fc='grey', alpha=0.7)
# ll_other3 = (xo3[j, 0]+rdiag*np.cos(a0+xo3[j, 2]),
#                  xo3[j, 1]+rdiag*np.sin(a0+xo3[j, 2]))
# show_other3 = patches.Rectangle(ll_other3, width, height,
#                               angle=xo3[j, 2]*180/math.pi,
#                               fill = True, ec='grey', fc='grey', alpha=0.7)
# ax.add_patch(show_ego)
# ax.add_patch(rshow_ego)
# ax.add_patch(show_other1)
# ax.add_patch(show_other2)
# ax.add_patch(show_other3)

plt.savefig(dirpath+'/traj_closedloop_ca_multi.pdf')

plt.show()
