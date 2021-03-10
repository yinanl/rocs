import numpy as np
from scipy.linalg import expm, inv
from scipy.optimize import fsolve


def dcdc(t, x, u):
    xc = 70
    xl = 3
    rc = 0.005
    rl = 0.05
    r0 = 1
    vs = 1
    A1 = np.array([[-rl/xl, 0], [0, -1/(xc*(rc+r0))]])
    b1 = np.array([[vs/xl], [0]])
    A2 = np.array([[(-1/xl)*(rl + r0*rc/(r0+rc)), (-1/xl)*(r0/(r0+rc))],
                   [(1/xc)*(r0/(r0+rc)), (-1/xc)*(1/(r0+rc))]])
    b2 = b1
    I = np.eye(2)
    if(u == 1):
        return np.matmul(expm(A1*t), x) \
            + np.matmul(inv(A1), np.matmul((expm(A1*t)-I), b1))
    elif(u == 2):
        return np.matmul(expm(A2*t), x) \
            + np.matmul(inv(A2), np.matmul((expm(A2*t)-I), b2))
    else:
        raise ValueError('dcdc: Wrong input value.')


def car(t, x, u):
    alpha = np.arctan(np.tan(u[1])/2.)
    return np.array([u[0]*np.cos(alpha+x[2])/np.cos(alpha),
                     u[0]*np.sin(alpha+x[2])/np.cos(alpha),
                     u[0]*np.tan(u[1])])


def car2(t, x, u):
    return np.array([u[0]*np.cos(x[2]),
                     u[0]*np.sin(x[2]),
                     u[1]])


def MG(t, x, u):
    a = 1./3.5
    B = 2.0
    H = 0.18
    W = 0.25
    lc = 8.0
    cx = 1.0/lc
    cy = 1.0/(4*lc*B*B)
    aH = a+H
    H2 = H/(2.0*W*W*W)
    W2 = 3*W*W

    dx = cx * (aH+H2*(x[0]-W)*(W2-(x[0]-W)*(x[0]-W)) - x[1]) + u[0]
    dy = cy * (x[0] - u[1]*np.sqrt(x[1]))
    return np.array([dx, dy])


def MG_scaled(t, x, u):
    a = 1./3.5
    B = 2.0
    H = 0.18
    W = 0.25
    lc = 8.0
    cx = 1.0/lc
    cy = 1.0/(4*lc*B*B)
    aH = a+H
    H2 = H/(2.0*W*W*W)
    W2 = 3*W*W

    s = 10.0
    dx = cx * (aH+H2*(x[0]-W)*(W2-(x[0]-W)*(x[0]-W)) - x[1]/s) + u[0]
    dy = s*cy * (x[0] - u[1]*np.sqrt(x[1]/s))
    return np.array([dx, dy])


def MG_EP(psi_c0, H, W, mu):
    psi_c = lambda x: psi_c0+H*(1+1.5*(x/W-1)-0.5*(np.power((x/W-1),3)))-x**2/mu/mu
    xe = fsolve(psi_c, 0.5)
    ye = np.square(xe/mu)
    return (xe, ye)


def MG_plotEP(ax):
    H = 0.18
    W = 0.25
    xe = []
    ye = []
    mu_0 = 0.5
    for i in range(0, 20):
        mu2 = mu_0+i*0.01
        x, y = MG_EP(1.67*H, H, W, mu2)

        xe.append(x[0])
        ye.append(y[0])
    ax.plot(xe[0:12],ye[0:12],'--', label='Ustable Equilibrium Points')
    ax.plot(xe[11:20],ye[11:20], label='Stable Equilibrium Points')
    ax.legend(loc='lower right')


def scara_2dbint(t, x, u):
    return np.array([x[2], x[3], u[0], u[1]])


class scara:
    def __init__(self, tau):
        self.tau = tau
        self.m1 = 0.1
        self.m2 = 0.1
        self.l1 = 0.15
        self.l2 = 0.15
        self.r1 = 0.5*self.l1
        self.r2 = 0.5*self.l2
        self.I1 = 1.33e-5
        self.I2 = 1.33e-5

        self.z1 = self.I1+self.I2+self.m1*self.r1**2 \
            + self.m2*(self.l1**2+self.r2**2)
        self.z2 = self.m2*self.l1*self.r2
        self.z3 = self.I2+self.m2*self.r2**2

    def xy2theta(self, x, y):
        l3 = np.sqrt(x**2 + y**2)
        a = np.acos((l3**2+self.l1**2-self.l2**2) / (2*self.l1*l3))
        b = np.acos((l3**2+self.l2**2-self.l1**2) / (2*self.l2*l3))
        theta1 = [np.atan(y/x) + a, np.atan(y/x) - a]
        theta2 = [-(a+b), a+b]
        return np.array([theta1, theta2])

    def theta2xy(self, theta1, theta2):
        x = self.l1 * np.cos(theta1) + self.l2 * np.cos(theta1+theta2)
        y = self.l1 * np.sin(theta1) + self.l2 * np.sin(theta1+theta2)
        return np.array([x, y])

    def ode_full(self, t, x, u):
        s2 = np.sin(x[1])
        c2 = np.cos(x[1])
        D = self.z3*(self.z1-self.z3) - self.z2**2 * c2**2

        dw1dt = (self.z3*u[0] + self.z3*self.z2*s2*(2*x[2]+x[3])*x[3]
                 + (self.z3+self.z2*c2)*(self.z2*x[2]**2*s2-u[1])) / D
        dw2dt = ((self.z1+2*self.z2*c2)*(u[1]-self.z2*x[2]**2*s2)
                 - (self.z3+self.z2*c2)*(u[0]+self.z2*s2*(2*x[2]+x[3])*x[3]))/D
        return np.array([x[2], x[3], dw1dt, dw2dt])

    def ode_2dbint(self, t, x, u):
        return np.array([x[2], x[3], u[0], u[1]])

    def de_2dbint(self, t, x, u):
        return np.array([x[0] + t*x[2] + 0.5*t**2*u[0],
                         x[2] + t*u[0],
                         x[1] + t*x[3] + 0.5*t**2*u[1],
                         x[3] + t*u[1]])

    def compute_torque(self, u, w, theta):
        '''
        Calculate the torque for corresponding angles and angular velocities.
        u(2 x 1):
        '''
        ct2 = np.cos(theta[1])
        st2 = np.sin(theta[1])
        M = np.array([[self.z1+2*self.z2*ct2, self.z3+self.z2*ct2],
                      [self.z3+self.z2*ct2, self.z3]])
        C = np.array([[-self.z2*st2*w[1], -self.z2*st2*(w[0]+w[1])],
                      [self.z2*st2*w[0], 0]])
        return np.matmul(M, u)+np.matmul(C, w)
