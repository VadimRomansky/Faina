import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_synchrotron3():
    radiationFile1 = open('../css161010_99.dat','r')
    radiationFile2 = open('../css161010_162.dat','r')
    radiationFile3 = open('../css161010.dat','r')
    lines1 = radiationFile1.readlines()
    N1 = len(lines1)

    radiation1 = np.zeros([2,N1])
    for i in range(N1):
        s = lines1[i].split()
        radiation1[0,i] = float(s[0])
        radiation1[1,i] = float(s[1])

    lines2 = radiationFile2.readlines()
    N2 = len(lines2)

    radiation2 = np.zeros([2,N2])
    for i in range(N2):
        s = lines2[i].split()
        radiation2[0,i] = float(s[0])
        radiation2[1,i] = float(s[1])

    lines3 = radiationFile3.readlines()
    N3 = len(lines1)

    radiation3 = np.zeros([2,N3])
    for i in range(N3):
        s = lines3[i].split()
        radiation3[0,i] = float(s[0])
        radiation3[1,i] = float(s[3])

    cssFile1 = open('../examples_data/css_data/coppejans99.txt', 'r')
    cssFile2 = open('../examples_data/css_data/coppejans162.txt', 'r')
    cssFile3 = open('../examples_data/css_data/coppejans357.txt', 'r')
    lines1 = cssFile1.readlines()
    Ncss1= len(lines1)

    css1 = np.zeros([3, Ncss1])
    for i in range(Ncss1):
        s = lines1[i].split()
        css1[0, i] = float(s[4])
        css1[1, i] = float(s[6])
        css1[2, i] = float(s[7])

    lines2 = cssFile2.readlines()
    Ncss2 = len(lines2)

    css2 = np.zeros([3, Ncss2])
    for i in range(Ncss2):
        s = lines2[i].split()
        css2[0, i] = float(s[4])
        css2[1, i] = float(s[6])
        css2[2, i] = float(s[7])

    lines3 = cssFile3.readlines()
    Ncss3= len(lines3)

    css3 = np.zeros([3, Ncss3])
    for i in range(Ncss3):
        s = lines3[i].split()
        css3[0, i] = float(s[4])
        css3[1, i] = float(s[6])
        css3[2, i] = float(s[7])

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
    plt.plot(radiation1[0], radiation1[1], 'r', linewidth=4, label = '99 days')
    plt.plot(radiation2[0], radiation2[1], 'g', linewidth=4, label = '162 days')
    plt.plot(radiation3[0], radiation3[1], 'b', linewidth=4, label = '357 days')
    plt.errorbar(css1[0,:], css1[1,:], css1[2,:], ecolor = 'r', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)
    plt.errorbar(css2[0,:], css2[1,:], css2[2,:], ecolor = 'g', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)
    plt.errorbar(css3[0,:], css3[1,:], css3[2,:], ecolor = 'b', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)
    ax.legend(fontsize="25")
    #plt.show()
    plt.savefig('synchrotron3.png', bbox_inches='tight')