import numpy as np


def engineMG(t, x, u):
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
