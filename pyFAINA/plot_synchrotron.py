import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_synchrotron():
    radiationFile = open('../css161010.dat','r')
    lines = radiationFile.readlines()
    N = len(lines)

    radiation = np.zeros([2,N])
    for i in range(N):
        s = lines[i].split()
        radiation[0,i] = float(s[0])
        radiation[1,i] = float(s[1])

    cssFile = open('../examples_data/css_data/coppejans99.txt', 'r')
    #cssFile = open('../examples_data/css_data/coppejans357.txt', 'r')
    lines = cssFile.readlines()
    Ncss = len(lines)

    css = np.zeros([3, Ncss])
    for i in range(Ncss):
        s = lines[i].split()
        css[0, i] = float(s[4])
        css[1, i] = float(s[6])
        css[2, i] = float(s[7])

    #cssx1 = [1.5, 3.0, 6.1, 9.87]
    #cssy1 = [1.5, 4.3, 6.1, 4.2]
    #cssError1 = [0.1, 0.2, 0.3, 0.2]



    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    ax.set_xlabel(r'$\nu~GHz$', fontsize=40,fontweight='bold')
    #ax.set_xlabel(r'$\nu~Ггц$', fontsize=40,fontweight='bold')
    ax.set_ylabel(r'$F_{\nu}~mJy$', fontsize=40,fontweight='bold')
    #ax.set_ylabel(r'$F_{\nu}~мЯн$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    ax.set_xlim([0.1, 200])
    ax.set_ylim([0.02, 20])
    ax.set_xscale("log")
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])
    plt.plot(radiation[0], radiation[1], 'r', linewidth=4)
    plt.errorbar(css[0,:], css[1,:], css[2,:], ecolor = 'b', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)
    #plt.show()
    plt.savefig('synchrotron1.png', bbox_inches='tight')