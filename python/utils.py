import h5py
import numpy as np
import re
from functools import reduce
from matplotlib.collections import PolyCollection


def read_controller_itvl_from_h5(filename):
    tau = 0.01
    X = np.array([])
    U = np.array([])
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
    return tau, X, U, tag, pavings, ctlr


def read_controller_abst_from_h5(filename):
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
        goalset = f['G'][...]
        winids = f['WinSet'][...]
        ctlr = f['OptCtlr'][...]
        encode3 = f['encode3'][...]
        nts_ctrlr = f['nts_ctrlr'][...]
        q_prime = f['q_prime'][...]
    return tau, X, U, eta, xgrid, goalset, winids, \
        ctlr, encode3, nts_ctrlr, q_prime


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
    return n_dba, n_props, q0, acc, M.astype(int)


def polycoll_winset_itvl(winset):
    xys = winset[:, 0:4]
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
    return PolyCollection(verts, closed=True,
                          color='palegoldenrod', alpha=0.7)


def polycoll_winset_abst(winset, eta):
    xys = winset[:, 0:2]
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
    return PolyCollection(verts, closed=True,
                          color='palegoldenrod', alpha=0.7)


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
