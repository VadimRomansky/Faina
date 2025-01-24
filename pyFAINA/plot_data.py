from matplotlib import animation
from pylab import *
import numpy as np

def plot_data(filename, name, Ncol):
    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4

    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)

    file = open(filename, 'r')
    lines = file.readlines()
    Nraw = len(lines)

    data = np.zeros([Ncol+1, Nraw])

    for i in range(Nraw):
        for j in range(Ncol+1):
            data[j][i] = float(lines[i].split()[j])

    for i in range(Ncol):
        plt.plot(data[0],data[i+1])

    plt.savefig(name + '.png', bbox_inches='tight')