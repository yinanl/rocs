import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
import matplotlib as mpl
import re
import numpy as np
import random
import math
from scipy.integrate import solve_ivp
import sys
sys.path.insert(1, '/Users/yinan/Desktop/rocs/python/')
import utils
from odes import car2


dirpath = "/Users/yinan/Desktop/rocs/examples/collision-avoid"
tau = 0.3
width = 0.6
height = 0.5  # RobotCar size
dmin = 0.8 # minimum distance to the other robot
umin = np.array([-0.8, -0.8])
umax = np.array([0.8, 0.8])


fig = plt.figure()
ax = plt.axes()
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)

# # # Scenario 1
# tau = 0.3
# width = 0.4
# height = 0.3  # RobotCar size
# dmin = 1.2 # minimum distance to the other robot

# obstacles = [[1.6+5, 5.7+5, 4.0+5, 5.0+5],
#              [3+5, 5+5, 5+5, 8+5],
#              [4.3+5, 5.7+5, 1.8+5, 4.0+5],
#              [5.7+5, 8.5+5, 1.8+5, 2.5+5]]
# rects_obs = [patches.Rectangle((obstacles[i][0], obstacles[i][2]),
#                                obstacles[i][1]-obstacles[i][0],
#                                obstacles[i][3]-obstacles[i][2],
#                                linewidth=1,edgecolor='k',facecolor='k')
#               for i in range(len(obstacles))]

# goals = [[1+3, 2+3, 0.5+3, 2+3],
#          [0.5+3, 2.5+3, 7.5+5, 8.5+5],
#          [7.1+7, 9.1+7, 4.6+7, 6.4+7]]
# rects_gs = [patches.Rectangle((goals[i][0], goals[i][2]),
#                                goals[i][1]-goals[i][0],
#                                goals[i][3]-goals[i][2],
#                                linewidth=1,edgecolor='y',facecolor='y')
#             for i in range(len(goals))]


# for rect in rects_obs:
#     ax.add_patch(rect)

# for rect in rects_gs:
#     ax.add_patch(rect)

# ax.text((goals[0][0]+goals[0][1])/2.0, (goals[0][2]+goals[0][3])/2.0, "$a_1$")
# ax.text((goals[1][0]+goals[1][1])/2.0, (goals[1][2]+goals[1][3])/2.0, "$a_2$")
# ax.text((goals[2][0]+goals[2][1])/2.0, (goals[2][2]+goals[2][3])/2.0, "$a_3$")


# # Scenario 2
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
                               linewidth=1,edgecolor='y',fill=False)
            for i in range(len(goals))]
circ_gs = patches.Circle([8.0, 8.0], radius=0.8, linewidth=1,
                         edgecolor='y', fill=False)


for rect in rects_obs:
    ax.add_patch(rect)

for rect in rects_gs:
    ax.add_patch(rect)
ax.add_patch(circ_gs)

ax.text((goals[0][0]+goals[0][1])/2.0-0.5, (goals[0][2]+goals[0][3])/2.0, "amount")
ax.text((goals[1][0]+goals[1][1])/2.0-0.5, (goals[1][2]+goals[1][3])/2.0, "drop")
ax.text(7.8, 7.8, "count")


# # load nominal ego trajectory
egoname = "/traj_closedloop_nominal"
egofile = dirpath + egoname + ".txt"
x = []
u = []
with open(egofile, "r") as ego:
    lines = ego.readlines()
    for line in lines:
        line = line.strip()
        xu = re.split(';', line)
        x.append([float(e) for e in re.split(',', xu[0])])
        u.append([float(e) for e in re.split(',|\n', xu[1])])

x = np.array(x)
u = np.array(u)


# # generate obstacle trajectory
np.random.default_rng()
# j0 = 25+3
# j1 = 33 #125
# j2 = 60
# j3 = 90
# wp = np.array([[2.44, 10.8, 0.0],
#                [3.7, 10.8, 0.0],
#                [7.5, 14.8, 0.96896],
#                [12.3, 14.8, -0.6018],
#                [15.2, 6.0, -1.2525]])
# r1 = 4.61
# r2 = np.sqrt((wp[3,0]-wp[2,0])**2+(wp[3,1]-wp[2,1])**2)/2 * r1 / (wp[2,0]-wp[1,0])

N = 300 #148
t = 0
j = 0
xi = np.array([-5., -5., 0.])
uo = np.array([0., 0.])
uo1 = np.array([0.2, 0.])
uo2 = np.array([0.5, 0.])
uo3 = np.array([0.2, 0.])

tsim = []
# xsim = []
# usim = []
xsim1 = []
usim1 = []
xsim2 = []
usim2 = []
xsim3 = []
usim3 = []
for j in range(N):
    # # Scenario 1
    # if(j < j0):
    #     y = xi
    #     tt = tau
    # else:
    #     if(j == j0):
    #         xo = wp[0, :]
    #     if(j < j1):
    #         uo[0] = 0.7
    #         uo[1] = 0.0
    #     elif(j < j2):
    #         uo[0] = 0.6
    #         uo[1] = uo[0]/r1
    #     elif(j < j3):
    #         uo[0] = 0.7
    #         uo[1] = -uo[0]/r2
    #     else:
    #         uo[0] = 0.8
    #         uo[1] = 0
    #     sol = solve_ivp(car2, [0, tau], xo, method='RK45', args=(uo,))
    #     tt = sol.t[-1]
    #     y = sol.y[:, -1]

    # # # Scenario 2:
    # if(j < 5):
    #     if(j == 0):
    #         xo = np.array([0., 3., 0.])
    #     uo[0] = 0.33
    #     uo[1] = 0
    # elif(j < 10):
    #     uo[0] = 0.8
    #     uo[1] = 0.7
    # elif(j < 36):
    #     uo[0] = 0.7
    #     uo[1] = 0
    # elif(j < 60):
    #     uo[0] = 0.7
    #     uo[1] = -0.6
    # elif(j < 70):
    #     uo[0] = 0.7
    #     uo[1] = 0.6
    # else:
    #     uo[0] = 0.5
    #     uo[1] = 0.0

    # # # Scenario 2: test case 2
    # if(j < 27):
    #     if(j == 0):
    #         xo = np.array([4.5, 0.0, np.pi/2.0])
    #     uo[0] = 0.65
    #     uo[1] = 0
    # elif(j < 30):
    #     uo[0] = 0.4
    #     uo[1] = -0.3
    # else:
    #     uo[0] = 0.4
    #     uo[1] = 0.3

    # # # Scenario 3:
    # if(t<18.9):
    #     if(t == 0):
    #         xo = np.array([10.0, 6.3, np.pi])
    #     uo[0] = 0.3
    #     uo[1] = 0.0
    # elif(t < 21):
    #     uo[0] = 0.3
    #     uo[1] = 0.4
    # elif(t < 25):
    #     uo[0] = 0.3
    #     uo[1] = 0.0
    # elif(t < 28):
    #     uo[0] = -0.3
    #     uo[1] = 0.0
    # elif(t < 29):
    #     uo[0] = 0.3
    #     uo[1] = 0.0
    # elif(t < 30):
    #     uo[0] = -0.3
    #     uo[1] = 0.0
    # else:
    #     # rn = np.random.rand(2)
    #     # uo = umin+(umax-umin)*rn
    #     uo[0] = 0.3
    #     uo[1] = 0.0


    # sol = solve_ivp(car2, [0, tau], xo, method='RK45', args=(uo,))
    # tt = sol.t[-1]
    # y = sol.y[:, -1]

    # # convert angle into [-pi, pi]
    # if(y[2] > np.pi):
    #     y[2] -= 2*np.pi
    # if(y[2] < -np.pi):
    #     y[2] += 2*np.pi

    # tsim.append(t)
    # xsim.append(xo)
    # # print(uo)
    # usim.append(uo)
    # # print(usim)
    # xo = y
    # t += tt

    # # Scenario 4: multi-obstacles
    # if(t < 19):
    #     if(t == 0):
    #         xo1 = np.array([0.0, 7.0, -0.6])
    #         xo2 = np.array([4.7, 0.0, 1.9])
    #         xo3 = np.array([10.0, 6.4, np.pi])
    # elif(t < 22):
    #     uo1[0] = 0.3
    #     uo1[1] = -0.4
    # else:
    #     uo1[0] = 0.3
    #     uo1[1] = 0.0
    #     uo2[0] = 0.3
    #     uo2[1] = 0.4
    if(t < 14):
        if(t == 0):
            xo1 = np.array([3.0, 10.0, -1.3])
            xo2 = np.array([5.3, 0.0, 2.0])
            xo3 = np.array([10.0, 6.4, np.pi])
    elif(t < 31):
        uo2[0] = 0.5
        uo2[1] = 0.6
    elif(t < 35):
        uo3[0] = 0.2
        uo3[1] = 0.3
    else:
        uo3[0] = 0.2
        uo3[1] = 0.0

    sol1 = solve_ivp(car2, [0, tau], xo1, method='RK45', args=(uo1,))
    y1 = sol1.y[:, -1]
    if(y1[2] > np.pi):
        y1[2] -= 2*np.pi
    if(y1[2] < -np.pi):
        y1[2] += 2*np.pi
    xsim1.append(xo1)
    usim1.append(uo1)
    sol2 = solve_ivp(car2, [0, tau], xo2, method='RK45', args=(uo2,))
    y2 = sol2.y[:, -1]
    if(y2[2] > np.pi):
        y2[2] -= 2*np.pi
    if(y2[2] < -np.pi):
        y2[2] += 2*np.pi
    xsim2.append(xo2)
    usim2.append(uo2)
    sol3 = solve_ivp(car2, [0, tau], xo3, method='RK45', args=(uo3,))
    y3 = sol3.y[:, -1]
    if(y3[2] > np.pi):
        y3[2] -= 2*np.pi
    if(y3[2] < -np.pi):
        y3[2] += 2*np.pi
    xsim3.append(xo3)
    usim3.append(uo3)

    tt = sol1.t[-1]
    tsim.append(t)

    t += tt
    xo1 = y1
    xo2 = y2
    xo3 = y3



# xo = np.asarray(xsim)
# uo = np.asarray(usim)
xo1 = np.asarray(xsim1)
uo1 = np.asarray(usim1)
xo2 = np.asarray(xsim2)
uo2 = np.asarray(usim2)
xo3 = np.asarray(xsim3)
uo3 = np.asarray(usim3)


# # animation
rdiag = np.sqrt(width**2+height**2)/2.0
a0 = np.arctan2(height, width)+np.pi
# lowleft_ego = (x[0, 0]+rdiag*np.cos(a0+x[0, 2]),
#                x[0, 1]+rdiag*np.sin(a0+x[0, 2]))
# rec_ego = patches.Rectangle(lowleft_ego, width, height,
#                             angle=x[0, 2]*180/math.pi,
#                             fill = True, ec='b', fc='b', alpha=0)
# # ax.add_patch(mpl.transforms.Affine2D().rotate_deg_around(x[0, 0], x[0, 1],
# #                                                          x[0,2]*180/np.pi))
# r_ego = patches.Circle((x[0,0], x[0,1]), radius=dmin/2.0,
#                        color='r', fill=False, alpha=0)
# ax.add_patch(rec_ego)
# ax.add_patch(r_ego)

# lowleft_other = (xo[0, 0]+rdiag*np.cos(a0+xo[0, 2]),
#                  xo[0, 1]+rdiag*np.sin(a0+xo[0, 2]))
# rec_other = patches.Rectangle(lowleft_other, width, height,
#                               angle=xo[0, 2]*180/np.pi,
#                               fill = True, ec='grey', fc='grey', alpha=0)
# ax.add_patch(rec_other)
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
    # rec_ego.xy = (x[i, 0]+rdiag*np.cos(a0+x[i, 2]),
    #               x[i, 1]+rdiag*np.sin(a0+x[i, 2]))
    # rec_ego.angle = x[i, 2] * 180 / np.pi

    # r_ego.center = (x[i,0], x[i,1])
    # rec_ego.set_alpha(0.7)
    # r_ego.set_alpha(0.7)

    # rec_other.xy = (xo[i, 0]+rdiag*np.cos(a0+xo[i, 2]),
    #                 xo[i, 1]+rdiag*np.sin(a0+xo[i, 2]))
    # rec_other.angle = xo[i, 2]*180/np.pi
    # rec_other.set_alpha(0.7)
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

    return rec_other1, rec_other2, rec_other3, time_text


ani = animation.FuncAnimation(fig, animate, range(len(x)),
                              interval=tau*500, blit=True)


# # Plot CA trajectory
# ax.plot(x[:, 0], x[:, 1])
# ax.plot(x[0,0], x[0,1], marker='^', markerfacecolor='r')
# ax.plot(x[-1,0], x[-1,1], marker='v', markerfacecolor='g')

# ax.plot(xo[xo[:,0]>0, 0], xo[xo[:,0]>0, 1], 'r')
# ax.plot(xo[-1,0], xo[-1,1], marker='v', markerfacecolor='g')
ax.plot(xo1[xo1[:,0]>0, 0], xo1[xo1[:,0]>0, 1], 'r')
ax.plot(xo1[-1,0], xo1[-1,1], marker='v', markerfacecolor='g')
ax.plot(xo2[xo2[:,0]>0, 0], xo2[xo2[:,0]>0, 1], 'r')
ax.plot(xo2[-1,0], xo2[-1,1], marker='v', markerfacecolor='g')
ax.plot(xo3[xo3[:,0]>0, 0], xo3[xo3[:,0]>0, 1], 'r')
ax.plot(xo3[-1,0], xo3[-1,1], marker='v', markerfacecolor='g')

# ax.set_aspect('equal')

plt.show()
