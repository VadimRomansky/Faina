from matplotlib import animation
from pylab import *
import numpy as np

def plot_distributions():
    Efile0 = open('../examples_data/gamma0.5_theta0-90/Ee1.dat','r').readline().split()
    Ffile0= open('../examples_data/gamma0.5_theta0-90/Fs1.dat','r').readline().split()

    Efile1 = open('../examples_data/gamma0.5_theta0-90/Ee3.dat','r').readline().split()
    Ffile1 = open('../examples_data/gamma0.5_theta0-90/Fs3.dat','r').readline().split()

    Efile2 = open('../examples_data/gamma0.5_theta0-90/Ee8.dat','r').readline().split()
    Ffile2 = open('../examples_data/gamma0.5_theta0-90/Fs8.dat','r').readline().split()

    N = len(Efile0)

    E = np.zeros([3,N])
    F = np.zeros([3,N])
    dE = np.zeros([3,N])

    factor = 1.2

    for i in range(N):
        E[0,i] = 1+float(Efile0[i])*factor
        E[1,i] = 1+float(Efile1[i])*factor
        E[2,i] = 1+float(Efile2[i])*factor

        F[0, i] = float(Ffile0[i])
        F[1, i] = float(Ffile1[i])
        F[2, i] = float(Ffile2[i])

        if i == 0 :
            dE[0,i] = float(Efile0[i+1]) - float(Efile0[i])
            dE[1,i] = float(Efile1[i+1]) - float(Efile1[i])
            dE[2,i] = float(Efile2[i+1])- float(Efile2[i])
        else:
            dE[0,i] = float(Efile0[i]) - float(Efile0[i-1])
            dE[1,i] = float(Efile1[i]) - float(Efile1[i-1])
            dE[2,i] = float(Efile2[i]) - float(Efile2[i-1])

    norm = np.zeros([3])
    for i in range(N):
        for j in range(3):
            norm[j] = norm[j] + F[j][i]*dE[j][i]

    for i in range(N):
        for j in range(3):
            F[j][i] = F[j][i]/norm[j]

    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    #plt.subplots_adjust(bottom=0.12)
    #ax.tick_params(axis='both', which='major', labelsize=10)
    #ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.set_xlabel('$E/m_e c^2$', fontsize=40,fontweight='bold')
    ax.set_ylabel('$F(E)$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    #ax.set_xlim([1E15, 1E17])
    ax.set_xscale("log")
    extraticks=[1,100]
    plt.xticks(list(plt.xticks()[0]) + extraticks)
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])
    plt.plot(E[0], F[0], 'b', linewidth=4)
    plt.plot(E[1], F[1], 'r', linewidth=4)
    plt.plot(E[2], F[2], 'g', linewidth=4)
    ax.legend([r'$\theta = 10^{\circ}$', r'$\theta = 30^{\circ}$', r'$\theta = 80^{\circ}$'], fontsize="40")

    #plt.show()
    plt.savefig('electrons.png', bbox_inches='tight')