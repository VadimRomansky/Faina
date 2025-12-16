from matplotlib import animation
from pylab import *
import numpy as np

def plot_thick_regime():
    plt.rcParams.update({'font.size': 40})
    #plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 1

    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    filename1 = '../testbigSource.dat'
    filename2 = '../thickCompton.dat'

    file1 = open(filename1, 'r')
    lines1 = file1.readlines()
    Nraw1 = len(lines1)

    data1 = np.zeros([2, Nraw1])

    for i in range(Nraw1):
        for j in range(2):
            data1[j][i] = float(lines1[i].split()[j])

    plt.plot(data1[0],data1[1], linewidth=4,label='big source')

    file2 = open(filename2, 'r')
    lines2 = file2.readlines()
    Nraw2 = len(lines2)

    data2 = np.zeros([2, Nraw2])

    for i in range(Nraw2):
        for j in range(2):
            data2[j][i] = float(lines2[i].split()[j])

    plt.plot(data2[0],data2[1], linewidth=4, label = 'thick regime')

    ax.set_ylim([1E-5, 1])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('E [eV]', fontsize=30,fontweight='bold')
    ax.set_ylabel('F(E) [erg/ cm^2 s]', fontsize=30,fontweight='bold')
    ax.legend(fontsize="20")


    plt.savefig('thick_regime.png', bbox_inches='tight')