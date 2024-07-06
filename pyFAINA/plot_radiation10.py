import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_radiation10():
    radiationFile0 = open('../outputSynch0.dat','r')

    lines = radiationFile0.readlines()
    N0 = len(lines)
    radiation = np.zeros([2*10,N0])

    for j in range(10):
        radiationFile = open('../outputSynch' + str(j) + '.dat', 'r')
        lines = radiationFile.readlines()
        for i in range(N0):
            s = lines[i].split()
            radiation[2*j,i] = float(s[0])
            radiation[2*j+1,i] = float(s[1])


    cssFile0 = open('../examples_data/css_data/coppejans69.txt')
    cssFile1 = open('../examples_data/css_data/coppejans99.txt')
    cssFile2 = open('../examples_data/css_data/coppejans162.txt')
    cssFile3 = open('../examples_data/css_data/coppejans357.txt')


    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    ax.set_xlabel(r'$\nu~GHz$', fontsize=40,fontweight='bold')
    ax.set_ylabel(r'$F_{\nu}~mJy$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    #ax.set_xlim([0.1, 100])
    #ax.set_ylim([0.02, 20])
    ax.set_xscale("log")
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])
    colors = ['r','indianred','orangered','salmon','maroon','b','aqua','navy','blueviolet','skyblue']
    for i in range(10):
        plt.plot(radiation[2*i,:], radiation[2*i+1,:], colors[i], linewidth = 4, label = str(i))
    ax.legend(fontsize="25")
    #plt.show()
    plt.savefig('radiation10.png', bbox_inches='tight')