from matplotlib import animation
from pylab import *
import numpy as np

def plot_data4(filename1, filename2, filename3, filename4, name, Ncol, label1 = "", label2 = "", label3 = "", label4 = "", xlabel = "", ylabel = ""):
    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4

    factor = 2E-13

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


    file4 = open(filename4, 'r')
    lines4 = file4.readlines()
    Nraw4 = len(lines4)

    data4 = np.zeros([Ncol + 1, Nraw4])

    for i in range(Nraw4):
        for j in range(Ncol + 1):
            data4[j][i] = float(lines4[i].split()[j])

    for i in range(Ncol):
        plt.plot(data1[0],data1[i+1]*factor, linewidth=4, label = label1)

    for i in range(Ncol):
        plt.plot(data2[0],data2[i+1]*factor, linewidth=4, label = label2)

    for i in range(Ncol):
        plt.plot(data3[0],data3[i+1]*3*factor, linewidth=4, label = label3)

    for i in range(Ncol):
        plt.plot(data4[0],data4[i+1]*factor, linewidth=4, label = label4)

    Ne = 100
    E1 = np.linspace(500, 10000, Ne)
    F1 = np.zeros([Ne])
    A = 8.2E-9
    for i in range(Ne):
        currentE = E1[i]*1.6E-12
        F1[i] = A*currentE**(2-1.58)

    E2 = np.linspace(3000,30000, Ne)
    F2 = np.zeros([Ne])
    B = 4.6E-12
    for i in range(Ne):
        currentE = E2[i]*1.6E-12
        F2[i] = B*currentE**(2-1.999)

    plt.plot(E1, F1, linewidth = 4, label = 'XMM')
    plt.plot(E2, F2, linewidth = 4, label = 'NuSTAR')

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel, fontsize=40, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=40, fontweight='bold')
    ax.legend(fontsize="20")

    ax.set_xlim([1, 1E16])
    ax.set_ylim([1E-13, 2E-10])


    plt.savefig(name + '.png', bbox_inches='tight')