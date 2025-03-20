from matplotlib import animation
from pylab import *
import numpy as np

def plot_diffusionConvection(filename, name, Ncol, xscale = 'linear', yscale = 'linear'):
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

    #for i in range(Ncol):
    #    plt.plot(data[0],data[i+1], linewidth=4)

    plt.plot(data[0], data[1], 'r', linewidth = 4, label = 'D d2f/dx2 E = 30 TeV')
    plt.plot(data[0], data[2], 'r', linewidth = 4, linestyle = 'dashed', label = 'u df/dx E = 30 TeV')

    plt.plot(data[0], data[3], 'b', linewidth=4, label='D d2f/dx2 E = 100 TeV')
    plt.plot(data[0], data[4], 'b', linewidth=4, linestyle='dashed', label='u df/dx E = 100 TeV')

    plt.plot(data[0], data[5], 'g', linewidth=4, label='D d2f/dx2 E = 300 TeV')
    plt.plot(data[0], data[6], 'g', linewidth=4, linestyle='dashed', label='u df/dx E = 300 TeV')

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    ax.set_xlim([-1E20, -1E10])
    #ax.set_ylim([1E-26, 1E-19])

    ax.legend(fontsize="20")

    plt.savefig(name + '.png', bbox_inches='tight')