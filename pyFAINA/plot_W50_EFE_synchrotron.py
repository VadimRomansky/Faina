import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_W50_EFE_synchrotron(filename, name, factor = 1.0):
    radiationFile = open(filename,'r')
    lines = radiationFile.readlines()
    N = len(lines)

    radiation = np.zeros([4,N])
    for i in range(N):
        s = lines[i].split()
        radiation[0,i] = float(s[0])
        radiation[1,i] = factor*float(s[1])
        radiation[2,i] = factor*float(s[2])
        radiation[3,i] = factor*float(s[3])

    outFile = open('Wsynchandcompt.dat', 'w')
    for i in range(N):
        print(radiation[0,i], radiation[1,i], radiation[2,i], radiation[3,i], sep = ' ', file = outFile)

    outFile.close()


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
    ax.set_xlim([1E2, 1E6])
    ax.set_ylim([0.5E-13, 0.5E-11])
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()
    #plt.xticks([1E11,2E11, 3E11, 4E11,5E11, 6E11, 7E11, 8E11, 9E11,1E12,2E12, 3E12, 4E12,5E12, 6E12, 7E12, 8E12, 9E12,1E13,2E13, 3E13, 4E13,5E13,6E13,7E13,8E13,9E13,1E14,2E14,3E14,4E14,5E14,6E14,7E14,8E14,9E14])
    plt.xticks([1E2, 2E2, 3E2, 4E2, 5E2, 6E2, 7E2, 8E2, 9E2, 1E3, 2E3, 3E3, 4E3, 5E3, 6E3, 7E3, 8E3, 9E3, 1E4, 2E4, 3E4, 4E4, 5E4, 6E4, 7E4, 8E4, 9E4, 1E5, 2E5, 3E5, 4E5, 6E5, 7E5, 8E5, 9E5, 1E6])
    #plt.yticks([6E-11, 7E-11, 8E-11, 9E-11, 1E-10, 2E-10])

    plt.plot(radiation[0], radiation[2], 'r', linewidth=4, label = 'head')
    plt.plot(radiation[0], radiation[3], 'b', linewidth=4, label = 'cone')

    Ne = 10
    E1 = np.linspace(500, 10000, Ne)
    F1h = np.zeros([Ne])
    F1c = np.zeros([Ne])
    F1d = np.zeros([Ne])
    Ah = 1.985E-9
    Ac = 5.2E-10
    Ad = 5.89E-11
    for i in range(Ne):
        currentE = E1[i]*1.6E-12
        F1h[i] = Ah*currentE**(2-1.58)
        F1c[i] = Ac*currentE**(2-1.65)
        F1d[i] = Ad*currentE**(2-1.75)

    E2 = np.linspace(3000,30000, Ne)
    F2h = np.zeros([Ne])
    F2c = np.zeros([Ne])
    F2d = np.zeros([Ne])
    Bh = 1.125E-9
    Bc = 6.86E-13
    Bd = 4.34E-13
    for i in range(Ne):
        currentE = E2[i]*1.6E-12
        F2h[i] = Bh*currentE**(2-1.6)
        F2c[i] = Bc*currentE**(2.0 - 2.0)
        F2d[i] = Bd*currentE**(2.0 - 2.0)

    plt.plot(E1, F1h, linewidth = 4, label = 'head XMM', color = 'm', linestyle = 'dashed')
    plt.plot(E2, F2h, linewidth = 4, label = 'head NuSTAR', color = 'm', linestyle = 'dotted')

    plt.plot(E1, F1c, linewidth = 4, label = 'cone XMM', color = 'c', linestyle = 'dashed')
    plt.plot(E2, F2c, linewidth = 4, label = 'cone NuSTAR', color = 'c', linestyle = 'dotted')

    ax.legend(fontsize = "20")
    #plt.plot(radiation[0], radiation[2], 'b', linewidth=4)
    #plt.errorbar(cssx1, cssy1, cssError1, ecolor = 'b', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)
    #ax.legend([r'BremsstrahlungThermalEvaluator', r'BremsstrahlungEbaluator'], fontsize="20")
    #plt.show()
    plt.savefig(name + '.png', bbox_inches='tight')