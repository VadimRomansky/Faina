from matplotlib import animation
from pylab import *
import numpy as np

def plot_W50_compton3(filename1, filename2, filename3, name, Ncol, label1 = "", label2 = "", label3 = "", xlabel = "", ylabel = ""):
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

    me = 0.91E-27
    c = 2.998E10
    me_c2 = me*c*c

    lhaasoFile = open("../examples_data/W50/LHAASO.dat", 'r')
    lhaasoLines = lhaasoFile.readlines()
    lhaasoN = len(lhaasoLines)

    lhaaso = np.zeros([4, lhaasoN])
    lhaasoLimits = np.zeros([lhaasoN])
    for i in range(lhaasoN):
        s = lhaasoLines[i].split()
        lhaaso[0, i] = float(s[0])
        lhaaso[1, i] = float(s[1])
        lhaaso[2, i] = float(s[2]) - lhaaso[1, i]
        lhaaso[3, i] = lhaaso[1, i] - float(s[3])
        lhaasoLimits[i] = False
    lhaasoLimits[lhaasoN - 1] = True
    lhaaso[2, lhaasoN - 1] = 0.5 * lhaaso[1, lhaasoN - 1]

    hessFile = open("../examples_data/W50/HESS.dat", 'r')
    hessLines = hessFile.readlines()
    hessN = len(hessLines)

    hess = np.zeros([6, hessN])
    hessLimits = np.zeros([hessN])
    for i in range(hessN):
        s = hessLines[i].split()
        hess[0, i] = float(s[0])
        hess[1, i] = float(s[1])
        hess[2, i] = float(s[2]) - hess[1, i]
        hess[3, i] = hess[1, i] - float(s[3])
        hess[4, i] = float(s[4]) - hess[0, i]
        hess[5, i] = hess[0, i] - float(s[5])
        hessLimits[i] = False
    hessLimits[hessN - 1] = True
    hess[3, i] = 0.5 * hess[3, i]

    for i in range(Ncol):
        plt.plot(data1[0],data1[i+1]*0.5E-11, linewidth=4, label = label1)

    for i in range(Ncol):
        plt.plot(data2[0],data2[i+1]*0.5E-11, linewidth=4, label = label2)

    for i in range(Ncol):
        plt.plot(data3[0],data3[i+1]*1E-4, linewidth=4, label = label3)

    plt.errorbar(lhaaso[0, :], lhaaso[1, :], yerr=[lhaaso[3, :], lhaaso[2, :]], uplims=lhaasoLimits, ecolor='r',
                 elinewidth=3, linewidth=0, capsize=5, capthick=3, label='LHAASO')
    plt.errorbar(hess[0, :], hess[1, :], yerr=[hess[3, :], hess[2, :]], xerr=[hess[5, :], hess[4, :]],
                 uplims=hessLimits, ecolor='purple', elinewidth=3, linewidth=0, capsize=5, capthick=3, label="H.E.S.S.")

    #ax.set_ylim([1E-14, 1E-11])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel, fontsize=40, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=40, fontweight='bold')
    ax.legend(fontsize="40")

    #.set_xlim([1E-2, 1E3])
    #ax.set_ylim([1E-1, 2E3])


    plt.savefig(name + '.png', bbox_inches='tight')