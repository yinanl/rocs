import scipy.io as sio
from scipy.integrate import solve_ivp
import h5py
import numpy as np
import re
from functools import reduce
from matplotlib.collections import PolyCollection
from collections import namedtuple


CtlrAbst = namedtuple('CtlrAbst', ['ugrid', 'xgrid', 'ctlr',
                                   'encode3', 'nts_ctrlr', 'q_prime'])
CtlrItvl = namedtuple('CtlrItvl', ['ugrid', 'pavings', 'tag', 'ctlr'])
CtlrSpec = namedtuple('CtlrSpec', ['n_dba', 'n_props',
                                   'q0', 'acc', 'q_prime'])


def read_controller_itvl_from_h5(filename):
    '''
    Make sure that the file contains at least:
    ts, X, U, tag, pavings, ctlr
    '''
    tau = 0.01
    X = np.array([])
    U = np.array([])
    G = np.array([])
    A = np.array([])
    tag = np.array([])
    pavings = np.array([])
    ctlr = np.array([])
    with h5py.File(filename, "r") as f:
        tau = f['ts'][...][0]
        X = f['X'][...]
        U = f['U'][...]
        tag = f['tag'][...]
        pavings = f['pavings'][...]
        ctlr = f['ctlr'][...]
        if("/G" in f):
            G = f['G'][...]
        if("/xobs" in f):
            A = f['xobs'][...]
    return tau, X, U, G, A, pavings, tag, ctlr


def read_controller_itvl_from_mat(filename):
    '''
    Make sure that the file contains at least:
    ts, X, U, tag, pavings, ctlr
    '''
    mat = sio.loadmat(filename)
    if("G" in mat):
        G = mat['G']
    else:
        G = np.array([])
    return mat['ts'][0, 0], mat['X'], mat['U'], G, \
        mat['tag'].squeeze(), mat['pavings'], mat['ctlr']


def read_controller_abst_from_h5(filename):
    '''
    Make sure that the file contains at least:
    ts, X, U, eta, xgrid, WinSet, OptCtlr, encode3, nts_ctrlr, and q_prime
    '''
    tau = np.array([])
    eta = np.array([])
    X = np.array([])
    U = np.array([])
    xgrid = np.array([])
    goalset = np.array([])
    winids = np.array([])
    ctlr = np.array([])
    encode3 = np.array([])
    nts_ctrlr = np.array([])
    q_prime = np.array([])
    with h5py.File(filename, "r") as f:
        tau = f['ts'][...][0]
        X = f['X'][...]
        U = f['U'][...]
        eta = f['eta'][...]
        xgrid = f['xgrid'][...]
        if("/G" in f):
            goalset = f['G'][...]
        winids = f['WinSet'][...]
        ctlr = f['OptCtlr'][...]
        encode3 = f['encode3'][...]
        nts_ctrlr = f['nts_ctrlr'][...]
        q_prime = f['q_prime'][...]
    # return tau, X, U, eta, xgrid, goalset, winids, \
    #     ctlr, encode3, nts_ctrlr, q_prime
    return tau, X, eta, goalset, winids, \
        CtlrAbst(U, xgrid, ctlr, encode3, nts_ctrlr, q_prime)


def read_spec_from_txt(filename):
    print("\nReading specification file...")
    n_dba = 0
    n_props = 0
    q0 = 0
    acc = 0
    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            cells = re.split('=', line)
            if(cells[0] == 'Name'):
                print(line)
            if(cells[0] == 'AP'):
                print(line)
            if(cells[0] == 'nNodes'):
                n_dba = int(cells[1])
            if(cells[0] == 'nAP'):
                n_props = 2**int(cells[1])
            if(cells[0] == 'init'):
                q0 = int(cells[1])
            if(cells[0] == 'acc'):
                acc = int(cells[1])
    M = np.zeros(shape=(n_dba, n_props))
    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            cells = re.split('#|,', line)
            if(len(cells) == 3):
                M[int(cells[0]), int(cells[1])] = int(cells[2])
    return CtlrSpec(n_dba, n_props, q0, acc, M.astype(int))


def polycoll_interval_array(arr, closed, color, alpha):
    '''
    Compile a 2d interval array into a collection of polytopes.
    arr (N x 2d): an array of size N with d-dimension (2d columns, [lower, upper]).
    closed: True or False.
    color: color for the polytopes.
    alpha: [0, 1] transparency.
    '''
    xys = arr[:, 0:4]
    tmp = np.zeros((xys.shape[0], xys.shape[1]*2))
    tmp[:, 0] = xys[:, 0]
    tmp[:, 1] = xys[:, 2]
    tmp[:, 2] = xys[:, 1]
    tmp[:, 3] = xys[:, 2]
    tmp[:, 4] = xys[:, 1]
    tmp[:, 5] = xys[:, 3]
    tmp[:, 6] = xys[:, 0]
    tmp[:, 7] = xys[:, 3]
    verts = tmp.reshape(xys.shape[0], 4, 2)
    return PolyCollection(verts, closed=closed, color=color, alpha=alpha)


def polycoll_grid_array(arr, eta, closed, color, alpha):
    '''
    Compile an array of data points into a collection of polytopes.
    arr (N x d): an array of size N with d-dimension (d columns).
    eta: an array of grid size (length of d).
    closed: True or False.
    color: color for the polytopes.
    alpha: [0, 1] transparency.
    '''
    xys = arr[:, 0:2]
    tmp = np.zeros((xys.shape[0], 8))
    tmp[:, 0] = xys[:, 0] - eta[0]/2.
    tmp[:, 1] = xys[:, 1] - eta[1]/2.
    tmp[:, 2] = xys[:, 0] + eta[0]/2.
    tmp[:, 3] = xys[:, 1] - eta[1]/2.
    tmp[:, 4] = xys[:, 0] + eta[0]/2.
    tmp[:, 5] = xys[:, 1] + eta[1]/2.
    tmp[:, 6] = xys[:, 0] - eta[0]/2.
    tmp[:, 7] = xys[:, 1] + eta[1]/2.
    verts = tmp.reshape(xys.shape[0], 4, 2)
    return PolyCollection(verts, closed=closed, color=color, alpha=alpha)


def index_in_interval_array(x, arr):
    '''
    Get the index of x in an interval array:
    x (n x 1): the input n-dim data point
    arr (m x 2n): an array of the n-dim intervals ([lower, upper])

    Returns: the first index where x belongs.
    '''
    if(arr.shape[1] < x.shape[0]*2):
        print("index_in_array: wrong dimension in input arrays.\n")
        return -1
    cond_left = np.asarray([x[i] >= arr[:, 2*i] for i in range(x.shape[0])])
    cond_right = np.asarray([x[i] <= arr[:, 2*i+1] for i in range(x.shape[0])])
    test_array = np.concatenate((cond_left, cond_right), axis=0)
    ids = np.argwhere(test_array.all(axis=0))
    if(not ids.size):
        return ids
    else:
        return ids[0, 0]


def index_in_grid(x, arr):
    '''
    Get the index of x in a grid arr:
    x (n x 1): the input n-dim data point
    arr (m x n): the n-dim grid (with m uniform data points)

    Returns: the index of the grid closest to x.
    '''
    d = x-arr
    d2 = [d[:, i]**2 for i in range(x.size)]
    dsum = reduce(lambda x, y: x+y, d2)
    return np.argmin(dsum)


def is_inside(x, w):
    '''
    Test if x is inside a given area w:
    x (n x 1): the input n-dim data point
    w (n x 2): the given n-dim interval area
    '''
    return np.all(x > w[:, 0]) and np.all(x < w[:, 1])


def simulate_abstbased_dba_control(tau, Tsim, num_acc, x0, vf, dba, controller):
    '''
    Generate sequences of closed-loop trajectories and control values.
    tau: Simulation time step.
    Tsim: Simulation time horizon.
    num_acc: The number of times visiting accepting nodes.
    x0: The initial condition.
    vf: Vector field of the dynamical system (or ODEs)
    dba: The CtlrSpec tuple with DBA info.
    controller: The CtlrAbst tuple with abstraction-based controller info.
    '''
    rng = np.random.default_rng()
    t = 0
    x = x0
    q = dba.q0
    tsim = []
    xsim = []
    usim = []
    qsim = []
    nacc = 0
    while(t < Tsim or nacc < num_acc):
        if(q == dba.acc):
            nacc += 1

        i = index_in_grid(x, controller.xgrid)  # convert x to node id
        p5 = controller.nts_ctrlr[controller.encode3[i], :]
        p7 = controller.ctlr[p5[2]:p5[2]+p5[0], :]
        uset = np.argwhere(p7[:, 0] == q).squeeze()
        if(uset.size > 1):
            uid = rng.choice(uset, 1)  # randomly pick one
        else:
            uid = int(uset)
            u = controller.ugrid[p7[uid, 1], :].squeeze()

        # Integrate ode
        sol = solve_ivp(vf, [0, tau], x, method='RK45', args=(u,))
        tt = sol.t[-1]
        y = sol.y[:, -1]
        if(y[2] > np.pi):
            y[2] -= 2*np.pi
        if(y[2] < -np.pi):
            y[2] += 2*np.pi

        # Update DBA state
        q = controller.q_prime[p5[1]*dba.n_dba+q]  # p5[1] is the label of current x

        # Save trajectories
        tsim.append(t)
        xsim.append(x)
        usim.append(u)
        qsim.append(q)

        # Update state
        x = y
        t += tt

    return np.asarray(xsim), np.asarray(usim), \
        np.asarray(qsim), np.asarray(tsim)


def simulate_itvl_dba_control(tau, Tsim, num_acc, x0, vf, dba,
                              controller, labeling):
    rng = np.random.default_rng()
    t = 0
    x = x0
    q = dba.q0
    tsim = []
    xsim = []
    usim = []
    qsim = []
    nacc = 0
    while(t < Tsim or nacc < num_acc):
        if(q == dba.acc):
            nacc += 1

        x_id = index_in_interval_array(x, controller.pavings[q])
        if(x_id < 0):
            print("System state ")
            print(x)
            print(" is not inside the winning set.")
            break
        if(q < 0):
            print("System unsafe.")
            break

        if(any(controller.ctlr[q][x_id, :])):
            uset = np.argwhere(controller.ctlr[q][x_id, :]).squeeze()  # get the indices of valid input
        else:
            print("No valid control input.")
            break

        if(uset.size > 1):
            uid = rng.choice(uset, 1)  # randomly pick one
        else:
            uid = int(uset)
        u = controller.ugrid[uid, :].squeeze()

        # Integrate ode
        sol = solve_ivp(vf, [0, tau], x, method='RK45', args=(u,))
        tt = sol.t[-1]
        y = sol.y[:, -1]
        if(y[2] > np.pi):
            y[2] -= 2*np.pi
        if(y[2] < -np.pi):
            y[2] += 2*np.pi

        # Save trajectories
        tsim.append(t)
        xsim.append(x)
        usim.append(u)
        qsim.append(q)
        # Update state
        q = dba.q_prime[q, labeling(x)]
        x = y
        t += tt

    return np.asarray(xsim), np.asarray(usim), \
        np.asarray(qsim), np.asarray(tsim)
