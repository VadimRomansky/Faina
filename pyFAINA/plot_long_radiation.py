import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_long_radiation():
    h = 6.626E-27
    tokeV = 1.6E-9

    radiationFile = open('../wideRangeSynch.dat','r')

    lines = radiationFile.readlines()
    N = len(lines)

    radiation = np.zeros([2,N])
    integral = 0;
    integral1 = 0;
    for i in range(N):
        s = lines[i].split()
        radiation[0,i] = float(s[0])
        if ((radiation[0, i] > 0.3*tokeV) and (radiation[0, i] < 10*tokeV) and (i > 0)):
            integral = integral + float(s[1]) * (radiation[0, i] - radiation[0, i - 1]*tokeV)

        if ((radiation[0, i] > 100*tokeV) and (radiation[0, i] < 1000*tokeV) and (i > 0)):
            integral1 = integral1 + float(s[1]) * (radiation[0, i] - radiation[0, i - 1]*tokeV)
        radiation[1,i] = 150*150*float(s[1])*radiation[0,i]
        radiation[0,i] = radiation[0,i]/tokeV

    radiationFile2 = open('../output3.dat','r')

    lines = radiationFile2.readlines()
    N2 = len(lines)

    radiation2 = np.zeros([2,N2])
    for i in range(N2):
        s = lines[i].split()
        radiation2[0,i] = float(s[0])
        radiation2[1,i] = 150*150*float(s[1])*radiation2[0,i]
        radiation2[0,i] = radiation2[0,i]/tokeV

    print("0.3-10 keV flux = ", integral)
    print("0.1-1 MeV flux = ", integral1)

    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    ax.set_xlabel('$\\rm E,~\text{keV}$', fontsize=40,fontweight='bold')
    #ax.set_ylabel('$\\rm E\,F(E),~10^{-12}~\text{эрг}~\text{см}^{-2}~\text{с{^{-1}$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    ax.set_xlim([1E-1, 0.5E4])
    #ax.set_ylim([1E-1, 2])
    ax.set_ylim([2, 20])
    ax.set_xscale("log")
    #extraticks=[1E-6,1E-2,100]
    #plt.xticks(list(plt.xticks()[0]+extraticks))
    #ax.set_xticks([1E-8, 1E-4, 1, 10000])
    ax.set_yticks([2, 5, 10, 20])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #ax1.get_xaxis().get_major_formatter().labelOnlyBase = False
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    #plt.yticks(fontsize=30)
    #plt.xticks(fontsize=30)
    #ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])

    plt.plot(radiation[0], radiation[1]*1E12, 'r', linewidth=4, label = r'синхротронное излучение')
    plt.plot(radiation2[0], radiation2[1]*1E12, 'b', linewidth=4, label = r'обратное комптоновское рассеяние')
    #plt.errorbar(cssx1, cssy1, cssError1, ecolor = 'b', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)
    ax.legend(fontsize="25")
    #plt.show()
    plt.savefig('long_radiation2.png', bbox_inches='tight')