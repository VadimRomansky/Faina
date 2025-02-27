from matplotlib import animation
from pylab import *
import numpy as np

def plot_data6(filename1, filename2, filename3, filename4, filename5, filename6, name, Ncol, label1 = "", label2 = "", label3 = "", label4 = "", label5 = "", label6 = "", xlabel = "", ylabel = ""):
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

    file5 = open(filename5, 'r')
    lines5 = file5.readlines()
    Nraw5 = len(lines5)

    data5 = np.zeros([Ncol + 1, Nraw5])

    for i in range(Nraw5):
        for j in range(Ncol + 1):
            data5[j][i] = float(lines5[i].split()[j])

    file6 = open(filename6, 'r')
    lines6 = file6.readlines()
    Nraw6 = len(lines6)

    data6 = np.zeros([Ncol + 1, Nraw6])

    for i in range(Nraw6):
        for j in range(Ncol + 1):
            data6[j][i] = float(lines6[i].split()[j])

    lhaasoFile = open("../examples_data/W50/LHAASO.dat",'r')
    lhaasoLines = lhaasoFile.readlines()
    lhaasoN = len(lhaasoLines)

    lhaaso = np.zeros([4, lhaasoN])
    lhaasoLimits = np.zeros([lhaasoN])
    for i in range(lhaasoN):
        s = lhaasoLines[i].split()
        lhaaso[0, i] = float(s[0])
        lhaaso[1, i] = float(s[1])
        lhaaso[2, i] = float(s[2]) - lhaaso[1,i]
        lhaaso[3, i] = lhaaso[1, i] - float(s[3])
        lhaasoLimits[i] = False
    lhaasoLimits[lhaasoN-1] = True
    lhaaso[2, lhaasoN - 1] = 0.5*lhaaso[1, lhaasoN - 1]

    xmmFile = open("../examples_data/W50/xmm_east.dat", 'r')
    xmmLines = xmmFile.readlines()
    xmmN = len(xmmLines)

    xmm = np.zeros([6, xmmN])
    for i in range(xmmN):
        s = xmmLines[i].split()
        xmm[0, i] = float(s[0])
        xmm[1, i] = float(s[1])
        xmm[2, i] = float(s[2]) - xmm[1, i]
        xmm[3, i] = xmm[1,i] - float(s[3])
        xmm[4, i] = float(s[4]) - xmm[0,i]
        xmm[5, i] = xmm[0,i] - float(s[5])

    fermiFile = open("../examples_data/W50/Fermi.dat", 'r')
    fermiLines = fermiFile.readlines()
    fermiN = len(fermiLines)

    fermi = np.zeros([5, fermiN])
    fermiLimits = np.zeros([fermiN])
    for i in range(fermiN):
        s = fermiLines[i].split()
        fermi[0, i] = float(s[0])
        fermi[1, i] = float(s[1])
        fermi[2, i] = float(s[2]) - fermi[0, i]
        fermi[3, i] = fermi[0, i] - float(s[3])
        fermi[4, i] = 0.5*fermi[1, i]
        fermiLimits[i] = True

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
    hess[3,i] = 0.5*hess[3,i]

    for i in range(Ncol):
        plt.plot(data1[0],data1[i+1]*factor, linewidth=4, label = label1)

    for i in range(Ncol):
        plt.plot(data2[0],data2[i+1]*factor, linewidth=4, label = label2)

    for i in range(Ncol):
        plt.plot(data3[0],data3[i+1]*factor, linewidth=4, label = label3)

    for i in range(Ncol):
        plt.plot(data4[0],data4[i+1]*factor, linewidth=4, label = label4)

    for i in range(Ncol):
        plt.plot(data5[0],data5[i+1]*factor, linewidth=4, label = label5)

    for i in range(Ncol):
        plt.plot(data6[0],data6[i+1]*factor, linewidth=4, label = label6)

    plt.errorbar(lhaaso[0, :], lhaaso[1, :], yerr = [lhaaso[3, :], lhaaso[2, :]], uplims = lhaasoLimits, ecolor='b', elinewidth=3, linewidth=0, capsize=5, capthick=3, label = 'LHAASO')
    plt.errorbar(xmm[0,:], xmm[1, :], yerr = [xmm[3,:], xmm[2, :]], xerr = [xmm[5, :], xmm[4, :]], ecolor = 'g', elinewidth=3, linewidth=0, capsize=5, capthick=3, label = 'XMM')
    plt.errorbar(fermi[0,:], fermi[1,:], yerr = fermi[4, :], xerr = [fermi[3, :], fermi[2, :]], uplims=fermiLimits, ecolor = 'y', elinewidth=3, linewidth=0, capsize=5, capthick=3, label = 'Fermi-LAT')
    plt.errorbar(hess[0,:], hess[1, :], yerr = [hess[3, :], hess[2, :]], xerr = [hess[5, :], hess[4, :]], uplims = hessLimits, ecolor='purple', elinewidth=3, linewidth=0, capsize=5, capthick=3, label = "H.E.S.S.")

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

    #ax.set_xlim([1, 1E17])
    #ax.set_ylim([1E-13, 2E-9])


    plt.savefig(name + '.png', bbox_inches='tight')