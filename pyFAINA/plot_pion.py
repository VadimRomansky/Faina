import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_pion():
    radiationFile = open('../outputPionE.dat','r')
    lines = radiationFile.readlines()
    N = len(lines)

    radiation = np.zeros([3,N])
    toTeV = 1.6E-12*1E12
    for i in range(N):
        s = lines[i].split()
        radiation[0,i] = float(s[0])/toTeV
        radiation[1,i] = float(s[1])

    fermiFile = open('../examples_data/cygnus/fermi.txt')
    lines = fermiFile.readlines()
    N = len(lines)
    fermi = np.zeros([6,N])
    for i in range(N):
        s = lines[i].split()
        for j in range(6):
            fermi[j][i] = float(s[j])

    argoFile = open('../examples_data/cygnus/argo.txt')
    lines = argoFile.readlines()
    N = len(lines)
    argo = np.zeros([6,N])
    for i in range(N):
        s = lines[i].split()
        for j in range(6):
            argo[j][i] = float(s[j])

    hwkFile = open('../examples_data/cygnus/hwk.txt')
    lines = hwkFile.readlines()
    N = len(lines)
    hwk = np.zeros([6,N])
    for i in range(N):
        s = lines[i].split()
        for j in range(6):
            hwk[j][i] = float(s[j])



    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    #ax.set_xlabel(r'$\nu~GHz$', fontsize=40,fontweight='bold')
    ax.set_xlabel(r'$E~GeV$', fontsize=40,fontweight='bold')
    #ax.set_xlabel(r'$\nu~Ггц$', fontsize=40,fontweight='bold')
    #ax.set_ylabel(r'$F_{\nu}~mJy$', fontsize=40,fontweight='bold')
    ax.set_ylabel(r'$EF(E)~TeV~cm^{-2} s^{-1}$', fontsize=40,fontweight='bold')
    #ax.set_ylabel(r'$F_{\nu}~мЯн$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    #ax.set_xlim([1E3, 1E8])
    #ax.set_ylim([1E-13, 1E-9])
    ax.set_xscale("log")
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])
    plt.plot(radiation[0], radiation[1], 'r', linewidth=4, label='model')
    #plt.plot(fermi[0], fermi[1], 'b', linewidth = 4)
    plt.errorbar(fermi[0], fermi[1], xerr = [fermi[2], fermi[3]], yerr = [fermi[5],fermi[4]], elinewidth = 4, linewidth=0, capsize = 5, capthick = 4, color = 'g', marker = 'o', label='Fermi LAT')
    plt.errorbar(argo[0], argo[1], yerr = [argo[5], argo[4]], elinewidth=4, linewidth=0, capsize = 5, capthick = 4, color = 'b', marker = 'o', label='ARGO')
    plt.errorbar(hwk[0], hwk[1], yerr = [hwk[5], hwk[4]], elinewidth=4, linewidth=0, capsize = 5, capthick = 4, color = 'black',marker = 'o', label='HAWC')
    ax.legend()
    #plt.show()
    plt.savefig('pion.png', bbox_inches='tight')