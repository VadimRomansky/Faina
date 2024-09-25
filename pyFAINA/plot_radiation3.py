import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_radiation3():
    radiationFile = open('../css161010.dat','r')


    lines = radiationFile.readlines()
    N0 = len(lines)
    radiation = np.zeros([4,N0])
    for i in range(N0):
        s = lines[i].split()
        for j in range(4):
            radiation[j,i] = float(s[j])


    cssFile0 = open('../examples_data/css_data/coppejans99.txt')
    cssFile1 = open('../examples_data/css_data/coppejans162.txt')
    cssFile2 = open('../examples_data/css_data/coppejans357.txt')

    lines = cssFile0.readlines()
    N0 = len(lines)
    css0 = np.zeros([3, N0])
    for i in range(N0):
        s = lines[i].split()
        css0[0, i] = float(s[4])
        css0[1,i] = float(s[6])
        css0[2,i] = float(s[7])

    lines = cssFile1.readlines()
    N1 = len(lines)
    css1 = np.zeros([3, N1])
    for i in range(N1):
        s = lines[i].split()
        css1[0, i] = float(s[4])
        css1[1,i] = float(s[6])
        css1[2,i] = float(s[7])

    lines = cssFile2.readlines()
    N2 = len(lines)
    css2 = np.zeros([3, N2])
    for i in range(N2):
        s = lines[i].split()
        css2[0, i] = float(s[4])
        css2[1,i] = float(s[6])
        css2[2,i] = float(s[7])


    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    ax.set_xlabel(r'$\nu~GHz$', fontsize=40,fontweight='bold')
    ax.set_ylabel(r'$F_{\nu}~mJy$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    ax.set_xlim([0.1, 100])
    ax.set_ylim([0.02, 20])
    ax.set_xscale("log")
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])

    #plt.plot(radiation0[0], radiation0[1], 'r', linewidth=4, label = '69 days')
    plt.plot(radiation[0], radiation[1], 'r', linewidth=4, label = '99 days')
    plt.plot(radiation[0], radiation[2], 'g', linewidth=4, label = '162 days')
    plt.plot(radiation[0], radiation[3], 'b', linewidth=4, label = '357 days')
    plt.errorbar(css0[0], css0[1], css0[2], ecolor = 'r', elinewidth = 3, linewidth=0, capsize = 5, capthick = 3)
    plt.errorbar(css1[0], css1[1], css1[2], ecolor = 'g', elinewidth = 3, linewidth=0, capsize = 5, capthick = 3)
    plt.errorbar(css2[0], css2[1], css2[2], ecolor = 'b', elinewidth = 3, linewidth=0, capsize = 5, capthick = 3)
    ax.legend(fontsize="25")
    #plt.show()
    plt.savefig('radiation3.png', bbox_inches='tight')