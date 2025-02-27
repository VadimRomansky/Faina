from matplotlib import animation
from pylab import *
import numpy as np

def plot_data2(filename1, filename2, name, Ncol, label1 = "", label2 = "", xlabel = "", ylabel = ""):
    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4

    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)

    file1 = open(filename1, 'r')
    lines1 = file1.readlines()
    Nraw1 = len(lines1)

    data1 = np.zeros([Ncol+1, Nraw1])

    for i in range(Nraw1):
        for j in range(Ncol+1):
            data1[j][i] = float(lines1[i].split()[j])

    file2 = open(filename2, 'r')
    lines2 = file2.readlines()
    Nraw2 = len(lines2)

    data2 = np.zeros([Ncol + 1, Nraw2])

    for i in range(Nraw2):
        for j in range(Ncol + 1):
            data2[j][i] = float(lines2[i].split()[j])

    me = 0.91E-27
    c = 2.998E10
    me_c2 = me*c*c

    for i in range(Ncol):
        plt.plot(data1[0],data1[i+1], linewidth=4, label = label1)

    for i in range(Ncol):
        plt.plot(data2[0],data2[i+1], linewidth=4, label = label2)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel, fontsize=40, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=40, fontweight='bold')
    ax.legend(fontsize="40")

    #.set_xlim([1E-2, 1E3])
    #ax.set_ylim([1E-1, 2E3])


    plt.savefig(name + '.png', bbox_inches='tight')