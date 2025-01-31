from matplotlib import animation
from pylab import *
import numpy as np

def plot_data3(filename1, filename2, filename3, name, Ncol):
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

    file3 = open(filename3, 'r')
    lines3 = file3.readlines()
    Nraw3 = len(lines3)

    data3 = np.zeros([Ncol + 1, Nraw3])

    for i in range(Nraw3):
        for j in range(Ncol + 1):
            data3[j][i] = float(lines3[i].split()[j])

    for i in range(Ncol):
        plt.plot(data1[0],data1[i+1])

    for i in range(Ncol):
        plt.plot(data2[0],data2[i+1])

    for i in range(Ncol):
        plt.plot(data3[0],data3[i+1])

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlim([1E-5, 1E4])
    ax.set_ylim([1E-16, 5E2])


    plt.savefig(name + '.png', bbox_inches='tight')