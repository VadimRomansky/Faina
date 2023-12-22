from matplotlib import animation
from pylab import *
import numpy as np

def plot_distributions_protons():
    Efile0 = open('../examples_data/gamma1.5_theta0-90_protons/Ee1.dat','r').readline().split()
    Ffile0= open('../examples_data/gamma1.5_theta0-90_protons/Fs1.dat','r').readline().split()

    Efile1 = open('../examples_data/gamma1.5_theta0-90_protons/Ee3.dat','r').readline().split()
    Ffile1 = open('../examples_data/gamma1.5_theta0-90_protons/Fs3.dat','r').readline().split()

    Efile2 = open('../examples_data/gamma1.5_theta0-90_protons/Ee8.dat','r').readline().split()
    Ffile2 = open('../examples_data/gamma1.5_theta0-90_protons/Fs8.dat','r').readline().split()

    Efile0a = open('../examples_data/gamma0.5_theta0-90_protons/Ee1.dat', 'r').readline().split()
    Ffile0a = open('../examples_data/gamma0.5_theta0-90_protons/Fs1.dat', 'r').readline().split()

    Efile1a = open('../examples_data/gamma0.5_theta0-90_protons/Ee3.dat', 'r').readline().split()
    Ffile1a = open('../examples_data/gamma0.5_theta0-90_protons/Fs3.dat', 'r').readline().split()

    Efile2a = open('../examples_data/gamma0.5_theta0-90_protons/Ee8.dat', 'r').readline().split()
    Ffile2a = open('../examples_data/gamma0.5_theta0-90_protons/Fs8.dat', 'r').readline().split()

    N = len(Efile0)

    E = np.zeros([6,N])
    F = np.zeros([6,N])
    dE = np.zeros([6,N])
    P = np.zeros([6,N])

    mp = 1.67E-24
    me = 0.91E-27
    c = 3*10E10
    c2 = c*c

    m = mp

    factor = 1.0

    for i in range(N):
        E[0,i] = (float(Efile0[i])*factor + 1)*m*c2
        E[1,i] = (float(Efile1[i])*factor + 1)*m*c2
        E[2,i] = (float(Efile2[i])*factor + 1)*m*c2
        E[3, i] = (float(Efile0a[i]) * factor + 1)*m*c2
        E[4, i] = (float(Efile1a[i]) * factor + 1)*m*c2
        E[5, i] = (float(Efile2a[i]) * factor + 1)*m*c2

        F[0, i] = float(Ffile0[i])
        F[1, i] = float(Ffile1[i])
        F[2, i] = float(Ffile2[i])
        F[3, i] = float(Ffile0a[i])
        F[4, i] = float(Ffile1a[i])
        F[5, i] = float(Ffile2a[i])

    for i in range(N):
        for j in range(6):
            if i == 0 :
                dE[j,i] = E[j][1] - E[j][0]
            else:
                dE[j,i] = E[j][i] - E[j][i-1]

            P[j,i] = sqrt(E[j,i]*E[j,i] - m*m*c2*c2)/c

    norm = np.zeros([6])
    for i in range(N):
        for j in range(6):
            norm[j] = norm[j] + F[j][i]*dE[j][i]

    for i in range(N):
        for j in range(6):
            F[j][i] = F[j][i]/norm[j]

    #for i in range(N):
    #    for j in range(6):
    #        F[j][i] = F[j][i]*E[j,i]*E[j,i]

    for i in range(N):
        for j in range(6):
            F[j][i] = (F[j][i]*P[j][i]*P[j][i]*P[j][i]*c2/E[j][i])/(m*c)
            P[j][i] = P[j][i]/(m*c)

    norm = np.zeros([6])
    for i in range(N):
        for j in range(6):
            dP = 0
            if i == 0:
                dP = P[j][1] - P[j][0]
            else:
                dP = P[j][i] - P[j][i - 1]
            norm[j] = norm[j] + (F[j][i] / (P[j][i] * P[j][i])) * dP

    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    #plt.subplots_adjust(bottom=0.12)
    #ax.tick_params(axis='both', which='major', labelsize=10)
    #ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.set_xlabel('$p/m_p c$', fontsize=40,fontweight='bold')
    ax.set_ylabel('$f(p/m_p c)  (p/m_p c)^4$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    ax.set_xscale("log")
    #extraticks=[1,100]
    plt.xticks(list(plt.xticks()[0]))
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])
    plt.plot(P[0], F[0], 'b', linewidth=4)
    plt.plot(P[3], F[3], 'b', linewidth=4, linestyle='dashed')
    plt.plot(P[1], F[1], 'r', linewidth=4)
    plt.plot(P[4], F[4], 'r', linewidth=4, linestyle='dashed')
    plt.plot(P[2], F[2], 'g', linewidth=4)
    plt.plot(P[5], F[5], 'g', linewidth=4, linestyle='dashed')
    ax.legend([r'$\theta = 10^{\circ}, v = 0.75c$',r'$\theta = 10^{\circ}, v = 0.5c$',
               r'$\theta = 30^{\circ}, v = 0.75c$',r'$\theta = 30^{\circ}, v = 0.5c$',
               r'$\theta = 80^{\circ}, v = 0.75c$',r'$\theta = 80^{\circ}, v = 0.5c$'], fontsize="25")
    ax.set_xlim(xmin=0.3, xmax=2E1)
    ax.set_ylim(ymin=1E-2, ymax=1)
    #ax.set_xbound(lower=0.5, upper=2E3)
    #plt.show()
    plt.savefig('protons_momentum.png', bbox_inches='tight')