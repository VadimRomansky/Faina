import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_long_radiation():
    radiationFile = open('../outputSynch5.dat','r')
    lines = radiationFile.readlines()
    N = len(lines)

    h=6.626E-27
    tokeV = 1.6E-9

    radiation = np.zeros([2,N])
    for i in range(N):
        s = lines[i].split()
        radiation[0,i] = float(s[0])
        radiation[1,i] = float(s[1])*radiation[0,i]
        radiation[0,i] = radiation[0,i]/tokeV;



    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    ax.set_xlabel(r'$E~keV$', fontsize=40,fontweight='bold')
    ax.set_ylabel(r'$E F(E)~erg~cm^{-2}~s^{-1}$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    #ax.set_xlim([0.5, 100])
    #ax.set_ylim([0.2, 20])
    ax.set_xscale("log")
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])
    plt.plot(radiation[0], radiation[1], 'r', linewidth=4)
    #plt.errorbar(cssx1, cssy1, cssError1, ecolor = 'b', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)

    #plt.show()
    plt.savefig('long_radiation.png', bbox_inches='tight')