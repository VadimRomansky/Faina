from matplotlib import animation
from pylab import *
import numpy as np

def plot_data(filename, name, Ncol, xscale = 'linear', yscale = 'linear'):
    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 1

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
        plt.plot(data[0],data[i+1], linewidth=4)

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    #ax.set_xlim([1, 10])
    #ax.set_ylim([1E-5, 7E-5])

    ax.legend([r'parallel', r'normal'], fontsize="30")

    plt.savefig(name + '.png', bbox_inches='tight')