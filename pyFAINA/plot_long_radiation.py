import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_long_radiation():
    radiationFile = open('../wideRangeSynch.dat','r')
    lines = radiationFile.readlines()
    N = len(lines)

    h=6.626E-27
    tokeV = 1.6E-9

    radiation = np.zeros([2,N])
    integral = 0;
    for i in range(N):
        s = lines[i].split()
        radiation[0,i] = float(s[0])
        if ((radiation[0, i] > 0.3*tokeV) and (radiation[0, i] < 10*tokeV) and (i > 0)):
            integral = integral + float(s[1]) * (radiation[0, i] - radiation[0, i - 1]*tokeV)
        radiation[1,i] = float(s[1])*radiation[0,i]/4
        radiation[0,i] = radiation[0,i]/tokeV

    print("0.3-10 keV flux = ", integral)

    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    ax.set_xlabel(r'$E~keV$', fontsize=40,fontweight='bold')
    ax.set_ylabel(r'$E F(E)~erg~cm^{-2}~s^{-1}$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    ax.set_xlim([1E-8, 0.5E4])
<<<<<<< HEAD
    ax.set_ylim([2E-16, 5E-15])
=======
    ax.set_ylim([1E-16, 2E-14])
>>>>>>> ff832d4f006fc109f26a686965a2ea6bd8b85f60
    ax.set_xscale("log")
    #extraticks=[1E-6,1E-2,100]
    #plt.xticks(list(plt.xticks()[0]+extraticks))
    ax.set_xticks([1E-8, 1E-6, 1E-4, 1E-2, 1, 100, 10000])
    ax.set_yticks([5E-16, 1E-15, 2E-15, 5E-15])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #ax1.get_xaxis().get_major_formatter().labelOnlyBase = False
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    plt.yticks(fontsize=30)
    plt.xticks(fontsize=30)
    ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])

    plt.plot(radiation[0], radiation[1], 'r', linewidth=4)
    #plt.errorbar(cssx1, cssy1, cssError1, ecolor = 'b', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)

    #plt.show()
    plt.savefig('long_radiation.png', bbox_inches='tight')