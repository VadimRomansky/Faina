import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_W50_EFE(filename, name):
    radiationFile = open(filename,'r')
    lines = radiationFile.readlines()
    N = len(lines)

    factor = 0.2E-13

    radiation = np.zeros([3,N])
    for i in range(N):
        s = lines[i].split()
        radiation[0,i] = float(s[0])
        radiation[1,i] = factor*float(s[1])
        #radiation[2,i] = float(s[2])

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



    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    #ax.set_xlabel(r'$\nu~GHz$', fontsize=40,fontweight='bold')
    ax.set_xlabel(r'$E~eV$', fontsize=40,fontweight='bold')
    #ax.set_xlabel(r'$\nu~Ггц$', fontsize=40,fontweight='bold')
    #ax.set_ylabel(r'$F_{\nu}~mJy$', fontsize=40,fontweight='bold')
    ax.set_ylabel(r'$EF(E)~erg~cm^{-2} s^{-1}$', fontsize=40,fontweight='bold')
    #ax.set_ylabel(r'$F_{\nu}~мЯн$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    ax.set_xscale("log")
    #ax.set_xlim([1E12, 1E16])
    #ax.set_xlim([1E11, 1E15])
    #ax.set_xlim([1E3, 5E4])
    ax.set_ylim([6E-15, 5E-10])
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()
    #plt.xticks([1E11,2E11, 3E11, 4E11,5E11, 6E11, 7E11, 8E11, 9E11,1E12,2E12, 3E12, 4E12,5E12, 6E12, 7E12, 8E12, 9E12,1E13,2E13, 3E13, 4E13,5E13,6E13,7E13,8E13,9E13,1E14,2E14,3E14,4E14,5E14,6E14,7E14,8E14,9E14])
    #plt.xticks([1E3, 2E3, 3E3, 4E3, 5E3, 6E3, 7E3, 8E3, 9E3, 1E4, 2E4, 3E4, 4E4, 5E4])
    #plt.yticks([6E-11, 7E-11, 8E-11, 9E-11, 1E-10, 2E-10])

    plt.plot(radiation[0], radiation[1], 'r', linewidth=4)
    #plt.plot(lhaaso[0], lhaaso[1],'b', linewidth=4)
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

    ax.legend(fontsize = "20")
    #plt.plot(radiation[0], radiation[2], 'b', linewidth=4)
    #plt.errorbar(cssx1, cssy1, cssError1, ecolor = 'b', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)
    #ax.legend([r'BremsstrahlungThermalEvaluator', r'BremsstrahlungEbaluator'], fontsize="20")
    #plt.show()
    plt.savefig(name + '.png', bbox_inches='tight')