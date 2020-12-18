import numpy as np


def car(t, x, u):
    alpha = np.arctan(np.tan(u[1])/2.)
    return np.array([u[0]*np.cos(alpha+x[2])/np.cos(alpha),
                     u[0]*np.sin(alpha+x[2])/np.cos(alpha),
                     u[0]*np.tan(u[1])])


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


class scara:
    def __init__(self, tau, m1, m2, l1, l2, r1, r2, I1, I2):
        self.tau = tau
        self.m1 = m1
        self.m2 = m2
        self.l1 = l1
        self.l2 = l2
        self.r1 = r1
        self.r2 = r2
        self.I1 = I1
        self.I2 = I2

        self.z1 = self.I1+self.I2+self.m1*self.r1**2+self.m2*(self.l1**2+self.r2**2)
        self.z2 = self.m2*self.l1*self.r1
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
        D = self.z3*(self.z1-self.z3) - self.z2*self.z2*c2**2

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


def scara_2dbint(t, x, u):
    return np.array([x[2], x[3], u[0], u[1]])
