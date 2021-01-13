import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
import re
import numpy as np
import random
import math
from scipy.integrate import solve_ivp


tau = 0.3
width = 0.4
height = 0.3  # RobotCar size
dmin = 1.2 # minimum distance to the other robot


fig = plt.figure()
ax = plt.axes()
ax.set_xlim(0, 10+10)
ax.set_ylim(0, 10+10)

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

# # load nominal ego trajectory
egoname = "sim_closedloop_nominal"
folder = "/Users/yinan/Desktop/rocs/examples/collision-avoid/"
egofile = folder + egoname + ".txt"
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
def car(t, x, u):
    return np.array([u[0]*np.cos(x[2]), u[0]*np.sin(x[2]), u[1]])


np.random.default_rng()

j0 = 25+3
j1 = 33 #125
j2 = 60
j3 = 90

wp = np.array([[2.44, 10.8, 0.0],
               [3.7, 10.8, 0.0],
               [7.5, 14.8, 0.96896],
               [12.3, 14.8, -0.6018],
               [15.2, 6.0, -1.2525]])
r1 = 4.61
r2 = np.sqrt((wp[3,0]-wp[2,0])**2+(wp[3,1]-wp[2,1])**2)/2 * r1 / (wp[2,0]-wp[1,0])

N = 148
t = 0
j = 0
xi = np.array([-5., -5., 0.])
xo = xi
uo = np.array([0., 0.])

xsim = []
tsim = []
usim = []
for j in range(N):
    if(j < j0):
        y = xi
        tt = tau
    else:
        # if(j == j1):
        #     xo = np.array([8.08, 3.67, 1.7])
        #     uo = np.array([0.3, 0.0])
        # else:
        #     uo[0] = np.random.uniform(-0.8, 0.8, 1)
        #     uo[1] = np.random.uniform(-0.8, 0.8, 1)

        if(j == j0):
            xo = wp[0, :]
        if(j < j1):
            uo[0] = 0.7
            uo[1] = 0.0
        elif(j < j2):
            uo[0] = 0.6
            uo[1] = uo[0]/r1
        elif(j < j3):
            uo[0] = 0.7
            uo[1] = -uo[0]/r2
        else:
            uo[0] = 0.8
            uo[1] = 0
        sol = solve_ivp(car, [0, tau], xo, method='RK45', args=(uo,))
        tt = sol.t[-1]
        y = sol.y[:, -1]
    # convert angle into [-pi, pi]
    if(y[2] > np.pi):
        y[2] -= 2*np.pi
    if(y[2] < -np.pi):
        y[2] += 2*np.pi

    tsim.append(t)
    xsim.append(xo)
    # print(uo)
    usim.append(uo)
    # print(usim)
    xo = y
    t += tt

xo = np.asarray(xsim)
uo = np.asarray(usim)


# # animation
lowleft_ego = (x[0, 0] + height/2*math.cos(x[0, 2]-math.pi/2),
               x[0, 1] + height/2*math.sin(x[0, 2]-math.pi/2))
rec_ego = patches.Rectangle(lowleft_ego, width, height,
                            angle=x[0, 2]*180/math.pi,
                            fill = True, ec='b', fc='b', alpha=0)
r_ego = patches.Circle((x[0,0], x[0,1]), radius=dmin,
                       color='r', fill=False, alpha=0)
ax.add_patch(rec_ego)
ax.add_patch(r_ego)

lowleft_other = (xo[0, 0] + height/2*math.cos(xo[0, 2]-math.pi/2),
                 xo[0, 1] + height/2*math.sin(xo[0, 2]-math.pi/2))
rec_other = patches.Rectangle(lowleft_other, width, height,
                              angle=xo[0, 2]*180/math.pi,
                              fill = True, ec='grey', fc='grey', alpha=0)
ax.add_patch(rec_other)

time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes)


def animate(i):
    time_text.set_text(time_template % (i * tau))
    rec_ego.xy = (x[i, 0] + height/2.*math.cos(x[i, 2]-math.pi/2.),
                  x[i, 1] + height/2.*math.sin(x[i, 2]-math.pi/2.))
    rec_ego.angle = x[i, 2] * 180 / np.pi
    r_ego.center = (x[i,0], x[i,1])
    rec_ego.set_alpha(0.7)
    r_ego.set_alpha(0.7)

    rec_other.xy = (xo[i, 0] + height/2*math.cos(xo[i, 2]-math.pi/2),
                    xo[i, 1] + height/2*math.sin(xo[i, 2]-math.pi/2))
    rec_other.angle = xo[i, 2] * 180 / math.pi
    rec_other.set_alpha(0.7)

    return rec_ego, rec_other, r_ego, time_text


ani = animation.FuncAnimation(fig, animate, range(len(x)),
                              interval=tau*500, blit=True)


# # Plot CA trajectory
ax.plot(x[:, 0], x[:, 1])
ax.plot(x[0,0], x[0,1], marker='^', markerfacecolor='r')
ax.plot(x[-1,0], x[-1,1], marker='v', markerfacecolor='g')
# ax.plot(x[76,0], x[76,1], marker='o', markerfacecolor='y')

ax.plot(xo[xo[:,0]>0, 0], xo[xo[:,0]>0, 1], 'r')
ax.plot(xo[-1,0], xo[-1,1], marker='v', markerfacecolor='g')

plt.show()
