import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_radiation_EFE(filename, name, xmin = None, xmax = None, ymin = None, ymax = None):
    radiationFile = open(filename,'r')
    lines = radiationFile.readlines()
    N = len(lines)

    factor = 2E-13

    radiation = np.zeros([3,N])
    for i in range(N):
        s = lines[i].split()
        radiation[0,i] = float(s[0])
        radiation[1,i] = float(s[1])*factor
        #radiation[2,i] = float(s[2])

    cssx1 = [1.5, 3.0, 6.1, 9.87]
    cssy1 = [1.5, 4.3, 6.1, 4.2]
    cssError1 = [0.1, 0.2, 0.3, 0.2]



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
    if((xmin != None) and (xmax != None)):
        ax.set_xlim([xmin, xmax])
    if((ymin != None) and (ymax != None)):
        ax.set_ylim([ymin, ymax])
    #ax.set_ylim([3E1, 8E1])
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()
    #plt.xticks([1E11,2E11,5E11,1E12,2E12,5E12,1E13,2E13,5E13,1E14,2E14,5E14])
    plt.xticks([1E3, 2E3, 3E3, 4E3, 5E3, 6E3, 7E3, 8E3, 9E3, 1E4, 2E4, 3E4, 4E4, 5E4, 6E4, 7E4, 8E4, 9E4, 1E5, 2E5])
    plt.yticks([2E-12, 3E-12, 4E-12, 5E-12, 6E-12, 7E-12, 8E-12, 9E-12, 1E-11, 2E-11, 3E-11, 4E-11, 5E-11, 6E-11, 7E-11, 8E-11, 9E-11, 1E-10, 2E-10])
    # plt.axis([0.0,1.0,0.0,1.0])
    plt.plot(radiation[0], radiation[1], 'r', linewidth=4)
    #plt.plot(radiation[0], radiation[2], 'b', linewidth=4)
    #plt.errorbar(cssx1, cssy1, cssError1, ecolor = 'b', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)
    #ax.legend([r'BremsstrahlungThermalEvaluator', r'BremsstrahlungEbaluator'], fontsize="20")
    #plt.show()
    plt.savefig(name + '.png', bbox_inches='tight')