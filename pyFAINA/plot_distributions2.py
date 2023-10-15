from matplotlib import animation
from pylab import *
import numpy as np

def plot_distributions2():
    Efile0 = open('../examples_data/gamma1.5_combined/Ee3.dat','r').readline().split()
    Ffile0= open('../examples_data/gamma1.5_combined/Fs3.dat','r').readline().split()

    Efile0a = open('../examples_data/gamma1.5_combined_protons/Ee3.dat', 'r').readline().split()
    Ffile0a = open('../examples_data/gamma1.5_combined_protons/Fs3.dat', 'r').readline().split()

    me = 9.1E-28
    mp = 1.6E-24
    c = 3E10

    N = len(Efile0)

    E = np.zeros([2,N])
    F = np.zeros([2,N])
    dE = np.zeros([2,N])


    for i in range(N):
        E[0,i] = float(Efile0[i])*me/mp
        E[1,i] = float(Efile0a[i])

        F[0, i] = float(Ffile0[i])
        F[1, i] = float(Ffile0a[i])

        for j in range(2):
            if i == 0:
                dE[j,i] = E[j,i+1]-E[j,i]
            else:
                dE[j,i] = E[j,i] - E[j,i-1]

    norm = np.zeros([2])
    for i in range(N):
        for j in range(2):
            norm[j] = norm[j] + F[j][i]*dE[j][i]

    Emax = (1.4E8)*me/mp
    for i in range(N):
        #F[0][i] = F[0][i]*exp(-E[0,i]/Emax)
        for j in range(2):
            F[j][i] = F[j][i]*E[j,i]*E[j,i]/norm[j]

    F3 = np.zeros([N])
    for i in range(N):
        F3[i] = F[0][i]*exp(-E[0,i]/Emax)

    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    #plt.subplots_adjust(bottom=0.12)
    #ax.tick_params(axis='both', which='major', labelsize=10)
    #ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.set_xlabel('$E_{kin}/m_p c^2$', fontsize=40,fontweight='bold')
    ax.set_ylabel('$F(E) E^2$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    ax.set_xscale("log")
    #extraticks=[1,100]
    #plt.xticks(list(plt.xticks()[0]))
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    #ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])
    plt.plot(E[0], F3, 'b', linewidth=4)
    plt.plot(E[1], F[1], 'r', linewidth=4)
    plt.plot(E[0], F[0], 'b', linewidth=4, linestyle='dashed')
    ax.legend([r'electrons',r'protons'], fontsize="25")
    #ax.set_xlim(xmin=0.5, xmax=2E3)
    ax.set_ylim(ymin=0.5E-6, ymax=2)
    #ax.set_xbound(lower=0.5, upper=2E3)
    #plt.show()
    plt.savefig('distributions.png', bbox_inches='tight')