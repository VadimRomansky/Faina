import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_W50_EFE(filename, name):
    radiationFile = open(filename,'r')
    lines = radiationFile.readlines()
    N = len(lines)

    factor = 2E-12

    radiation = np.zeros([3,N])
    for i in range(N):
        s = lines[i].split()
        radiation[0,i] = float(s[0])
        radiation[1,i] = factor*float(s[1])
        #radiation[2,i] = float(s[2])

    lhaasoFile = open("../examples_data/W50/LHAASO.dat",'r')
    lines1 = lhaasoFile.readlines()
    N1 = len(lines1)

    radiation1 = np.zeros([3, N1])
    for i in range(N1):
        s = lines1[i].split()
        radiation1[0, i] = float(s[0])
        radiation1[1, i] = float(s[1])
        radiation1[2, i] = float(s[2]) - radiation1[1,i]



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
    ax.set_xlim([1E11, 1E15])
    #ax.set_ylim([3E1, 8E1])
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()
    plt.xticks([1E11,2E11, 3E11, 4E11,5E11, 6E11, 7E11, 8E11, 9E11,1E12,2E12, 3E12, 4E12,5E12, 6E12, 7E12, 8E12, 9E12,1E13,2E13, 3E13, 4E13,5E13,6E13,7E13,8E13,9E13,1E14,2E14,3E14,4E14,5E14,6E14,7E14,8E14,9E14])
    #plt.xticks([1E3, 2E3, 3E3, 4E3, 5E3, 6E3, 7E3, 8E3, 9E3, 1E4, 2E4, 3E4, 4E4, 5E4])
    #plt.yticks([3E1, 4E1, 5E1, 6E1, 7E1, 8E1])
    # plt.axis([0.0,1.0,0.0,1.0])
    plt.plot(radiation[0], radiation[1], 'r', linewidth=4)
    #plt.plot(radiation1[0], radiation1[1],'b', linewidth=4)
    plt.errorbar(radiation1[0, :], radiation1[1, :], radiation1[2, :], ecolor='b', elinewidth=4, linewidth=0, capsize=5, capthick=4)
    #plt.plot(radiation[0], radiation[2], 'b', linewidth=4)
    #plt.errorbar(cssx1, cssy1, cssError1, ecolor = 'b', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)
    #ax.legend([r'BremsstrahlungThermalEvaluator', r'BremsstrahlungEbaluator'], fontsize="20")
    #plt.show()
    plt.savefig(name + '.png', bbox_inches='tight')