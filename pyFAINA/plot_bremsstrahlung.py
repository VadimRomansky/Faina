import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_bremsstrahlung():
    radiationFile = open('../outputBremNu.dat','r')
    lines = radiationFile.readlines()
    N = len(lines)
    Nrho = 5

    radiation = np.zeros([2 + 2*Nrho,N])
    for i in range(N):
        s = lines[i].split()
        radiation[0,i] = float(s[0])
        for j in range(Nrho):
            radiation[2*j + 1,i] = float(s[2*j+1])
            radiation[2*j + 2,i] = float(s[2*j+2])
        radiation[2*Nrho+1, i] = float(s[2*Nrho+1])



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
    #ax.set_xlim([1E3, 1E8])
    ax.set_ylim([1E4, 1E20])
    ax.set_xscale("log")
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])
    for j in range(Nrho):
        plt.plot(radiation[0], radiation[2*j+1], 'r', linewidth=4)
        plt.plot(radiation[0], radiation[2*j+2], 'b', linewidth=4)
    plt.plot(radiation[0], radiation[2*Nrho+1],'g', linewidth = 4)
    #plt.errorbar(cssx1, cssy1, cssError1, ecolor = 'b', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)
    #ax.legend([r'BremsstrahlungThermalEvaluator', r'BremsstrahlungEbaluator'], fontsize="20")
    #plt.show()
    plt.savefig('bremsstrahlung.png', bbox_inches='tight')