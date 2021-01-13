import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
import os
import re
import numpy as np
import math

fig = plt.figure()
ax = plt.axes()
ax.set_xlim(0,10+10)
ax.set_ylim(0,10+10)

# # Setup workspace
obstacles = [[1.6+5, 5.7+5, 4.0+5, 5.0+5],
             [3+5, 5+5, 5+5, 8+5],
             [4.3+5, 5.7+5, 1.8+5, 4.0+5],
             [5.7+5, 8.5+5, 1.8+5, 2.5+5]]
rects_obs = [patches.Rectangle((obstacles[i][0], obstacles[i][2]),
                               obstacles[i][1]-obstacles[i][0],
                               obstacles[i][3]-obstacles[i][2],
                               linewidth=1,edgecolor='k',facecolor='k')
              for i in range(len(obstacles))]

goals = [[1+3, 2+3, 0.5+3, 2+3],
         [0.5+3, 2.5+3, 7.5+5, 8.5+5],
         [7.1+7, 9.1+7, 4.6+7, 6.4+7]]
rects_gs = [patches.Rectangle((goals[i][0], goals[i][2]),
                               goals[i][1]-goals[i][0],
                               goals[i][3]-goals[i][2],
                               linewidth=1,edgecolor='y',facecolor='y')
            for i in range(len(goals))]


for rect in rects_obs:
    ax.add_patch(rect)

for rect in rects_gs:
    ax.add_patch(rect)

ax.text((goals[0][0]+goals[0][1])/2.0, (goals[0][2]+goals[0][3])/2.0, "$a_1$")
ax.text((goals[1][0]+goals[1][1])/2.0, (goals[1][2]+goals[1][3])/2.0, "$a_2$")
ax.text((goals[2][0]+goals[2][1])/2.0, (goals[2][2]+goals[2][3])/2.0, "$a_3$")

# # load simulated trajectories
egoname = "sim_closedloop_ca_3"
othername = "sim_other_robot_3"
folder = "/Users/yinan/Desktop/rocs/examples/collision-avoid/"
egonom = folder + "sim_closedloop_nominal.txt"
egofile = folder + egoname + ".txt"
otherfile = folder + othername + ".txt"
x = []
u = []
xo = []
with open(egofile, "r") as ego, open(otherfile, "r") as other:
    lines = ego.readlines()
    for line in lines:
        line = line.strip()
        xu = re.split(';', line)
        x.append([float(e) for e in re.split(',', xu[0])])
        u.append([float(e) for e in re.split(',|\n', xu[1])])
    olines = other.readlines()
    for line in olines:
        line = line.strip()
        xu = re.split(';', line)
        xo.append([float(e) for e in re.split(',', xu[0])])
        # uo.append([float(e) for e in re.split(',|\n', xu[1])])


x = np.array(x)
u = np.array(u)
xo = np.array(xo)
# print(x)
# print(u)
# print(xo)


# # animation
tau = 0.3
width = 0.4
height = 0.3  # RobotCar size
dmin = 1.2 # minimum distance to the other robot
lowleft_ego = (x[0, 0] + height/2*math.cos(x[0, 2]-math.pi/2),
               x[0, 1] + height/2*math.sin(x[0, 2]-math.pi/2))
rec_ego = patches.Rectangle(lowleft_ego, width, height,
                            angle=x[0, 2]*180/math.pi,
                            fill = True, ec='b', fc='b', alpha=0)
lowleft_other = (xo[0, 0] + height/2*math.cos(xo[0, 2]-math.pi/2),
                 xo[0, 1] + height/2*math.sin(xo[0, 2]-math.pi/2))
rec_other = patches.Rectangle(lowleft_other, width, height,
                              angle=xo[0, 2]*180/math.pi,
                              fill = True, ec='grey', fc='grey', alpha=0)
r_ego = patches.Circle((x[0,0], x[0,1]), radius=dmin,
                       color='r', fill=False, alpha=0)

ax.add_patch(rec_ego)
ax.add_patch(rec_other)
ax.add_patch(r_ego)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes)

def animate(i):
    time_text.set_text(time_template % (i * tau))
    rec_ego.xy = (x[i, 0] + height/2.*math.cos(x[i, 2]-math.pi/2.),
                  x[i, 1] + height/2.*math.sin(x[i, 2]-math.pi/2.))
    rec_ego.angle = x[i, 2] * 180 / np.pi
    rec_other.xy = (xo[i, 0] + height/2*math.cos(xo[i, 2]-math.pi/2),
                    xo[i, 1] + height/2*math.sin(xo[i, 2]-math.pi/2))
    rec_other.angle = xo[i, 2] * 180 / math.pi
    r_ego.center = (x[i,0], x[i,1])
    rec_ego.set_alpha(0.7)
    rec_other.set_alpha(0.7)
    r_ego.set_alpha(0.7)
    return rec_ego, rec_other, r_ego, time_text

ani = animation.FuncAnimation(fig, animate, range(len(x)),
                              interval=tau*500, blit=True)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
ani.save(folder+'case1.mp4', writer=writer)


# # Plot nominal trajectory
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

# xo_dis = xo[xo[:,0]>0, :]
# ax.plot(xo[xo[:,0]>0, 0], xo[xo[:,0]>0, 1], 'r')
# ax.plot(xo_dis[0, 0], xo_dis[0, 1], marker='^', markerfacecolor='r')
# ax.plot(xo[-1,0], xo[-1,1], marker='v', markerfacecolor='g')

plt.savefig(folder+egoname+'.pdf')

plt.show()
