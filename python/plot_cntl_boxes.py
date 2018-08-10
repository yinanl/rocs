#
# call this function from command line:
# python plot
# plot from .txt file
# data in the form of [x1l, x1u, x2l, x2u, tag, cntl]

import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.ion()

parser = argparse.ArgumentParser(description='Plot 2D interval boxes.')
parser.add_argument("fname", help="indicate the file to plot")
parser.add_argument("-r", "--rectangle", help="plot rectangles",
                    action="store_true", default=1)
parser.add_argument("-c", "--center", help="plot center points",
                    action="store_true", default=0)
args = parser.parse_args()

xl = []
xu = []
yl = []
yu = []

# filename = 'dcdc_inv.txt'
# filename = 'data_8p6_roa.txt'
# filename = 'data_8p6_core.txt'
# fname = sys.argv[1]

plt.figure()
plt.axes()
# read one line a time, and close the file when done
with open(args.fname, 'r') as f:
    # lines = f.readlines()  # this will read everything into memory
    for line in f:
        xy = line.split()
        controls = [int(e) for e in xy[5:]]
        if np.any(controls):  # plot the box if there exist a control value
            xl = float(xy[0])
            xu = float(xy[1])
            yl = float(xy[2])
            yu = float(xy[3])
            if args.center:
                plt.plot((xl + xu)/2.0, (yl + yu)/2.0, '.')
            if args.rectangle:
                plt.gca().add_patch(
                    patches.Rectangle(
                        (xl, yl), xu - xl, yu - yl, fc='y', ec='k'))

        # xl.append(float(box[0]))
        # xu.append(float(box[1]))
        # yl.append(float(box[2]))
        # yu.append(float(box[3]))
        # tag.append(int(box[4]))

plt.axis('tight')  # others: 'auto', 'tight', 'scaled', 'square'
# plt.axis([1.15, 1.55, 1.09, 1.17])
# plt.axis([-1.5, 1.5, -1.5, 1.5])

raw_input("Press enter to exit...")
