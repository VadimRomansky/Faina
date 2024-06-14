from matplotlib import animation
from pylab import *
import numpy as np

def plot_distributions9():
    Efile0 = open('../examples_data/gamma1.5_theta0-90/Ee0.dat','r').readline().split()
    Ffile0= open('../examples_data/gamma1.5_theta0-90/Fs0.dat','r').readline().split()

    Efile1 = open('../examples_data/gamma1.5_theta0-90/Ee1.dat','r').readline().split()
    Ffile1 = open('../examples_data/gamma1.5_theta0-90/Fs1.dat','r').readline().split()

    Efile2 = open('../examples_data/gamma1.5_theta0-90/Ee2.dat','r').readline().split()
    Ffile2 = open('../examples_data/gamma1.5_theta0-90/Fs2.dat','r').readline().split()

    Efile3 = open('../examples_data/gamma1.5_theta0-90/Ee3.dat', 'r').readline().split()
    Ffile3 = open('../examples_data/gamma1.5_theta0-90/Fs3.dat', 'r').readline().split()

    Efile4 = open('../examples_data/gamma1.5_theta0-90/Ee4.dat', 'r').readline().split()
    Ffile4 = open('../examples_data/gamma1.5_theta0-90/Fs4.dat', 'r').readline().split()

    Efile5 = open('../examples_data/gamma1.5_theta0-90/Ee5.dat', 'r').readline().split()
    Ffile5 = open('../examples_data/gamma1.5_theta0-90/Fs5.dat', 'r').readline().split()

    Efile6 = open('../examples_data/gamma1.5_theta0-90/Ee6.dat', 'r').readline().split()
    Ffile6 = open('../examples_data/gamma1.5_theta0-90/Fs6.dat', 'r').readline().split()

    Efile7 = open('../examples_data/gamma1.5_theta0-90/Ee7.dat', 'r').readline().split()
    Ffile7 = open('../examples_data/gamma1.5_theta0-90/Fs7.dat', 'r').readline().split()

    Efile8 = open('../examples_data/gamma1.5_theta0-90/Ee8.dat', 'r').readline().split()
    Ffile8 = open('../examples_data/gamma1.5_theta0-90/Fs8.dat', 'r').readline().split()

    Efile9 = open('../examples_data/gamma1.5_theta0-90/Ee9.dat', 'r').readline().split()
    Ffile9 = open('../examples_data/gamma1.5_theta0-90/Fs9.dat', 'r').readline().split()


    N = len(Efile0)

    E = np.zeros([10,N])
    F = np.zeros([10,N])
    dE = np.zeros([10,N])
    P = np.zeros([10,N])

    mp = 1.67E-24
    me = 0.91E-27
    c = 3*10E10
    c2 = c*c

    m = me

    factor = 1.2

    for i in range(N):
        E[0,i] = (float(Efile0[i])*factor + 1)*m*c2
        E[1,i] = (float(Efile1[i])*factor + 1)*m*c2
        E[2,i] = (float(Efile2[i])*factor + 1)*m*c2
        E[3, i] = (float(Efile3[i]) * factor + 1) * m * c2
        E[4, i] = (float(Efile4[i]) * factor + 1) * m * c2
        E[5, i] = (float(Efile5[i]) * factor + 1) * m * c2
        E[6, i] = (float(Efile6[i]) * factor + 1) * m * c2
        E[7, i] = (float(Efile7[i]) * factor + 1) * m * c2
        E[8, i] = (float(Efile8[i]) * factor + 1) * m * c2
        E[9, i] = (float(Efile9[i]) * factor + 1) * m * c2

        F[0, i] = float(Ffile0[i])
        F[1, i] = float(Ffile1[i])
        F[2, i] = float(Ffile2[i])
        F[3, i] = float(Ffile3[i])
        F[4, i] = float(Ffile4[i])
        F[5, i] = float(Ffile5[i])
        F[6, i] = float(Ffile6[i])
        F[7, i] = float(Ffile7[i])
        F[8, i] = float(Ffile8[i])
        F[9, i] = float(Ffile9[i])


    for i in range(N):
        for j in range(10):
            if i == 0 :
                dE[j,i] = E[j][1] - E[j][0]
            else:
                dE[j,i] = E[j][i] - E[j][i-1]

            P[j,i] = sqrt(E[j,i]*E[j,i] - m*m*c2*c2)/c

    norm = np.zeros([10])
    for i in range(N):
        for j in range(10):
            norm[j] = norm[j] + F[j][i]*dE[j][i]

    for i in range(N):
        for j in range(10):
            F[j][i] = F[j][i]/norm[j]

    #for i in range(N):
    #    for j in range(6):
    #        F[j][i] = F[j][i]*E[j,i]*E[j,i]

    for i in range(N):
        for j in range(10):
            F[j][i] = (F[j][i]*P[j][i]*P[j][i]*P[j][i]*c2/E[j][i])/(m*c)
            P[j][i] = P[j][i]/(m*c)

    norm = np.zeros([10])
    for i in range(N):
        for j in range(10):
            dP = 0
            if i == 0:
                dP = P[j][1] - P[j][0]
            else:
                dP = P[j][i] - P[j][i-1]
            ratio1 = dP/dE[j][i]
            ratio2 = E[j][i]/(P[j][i]*m*m*c2*c2)
            norm[j] = norm[j] + (F[j][i]/(P[j][i]*P[j][i])) * dP

    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    #plt.subplots_adjust(bottom=0.12)
    #ax.tick_params(axis='both', which='major', labelsize=10)
    #ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.set_xlabel('$p/m_e c$', fontsize=40,fontweight='bold')
    ax.set_ylabel('$f(p/m_e c)  (p/m_e c)^4$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    ax.set_xscale("log")
    #extraticks=[1,100]
    plt.xticks(list(plt.xticks()[0]))
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])
    plt.plot(P[0], F[0], 'b', linewidth=4)
    plt.plot(P[1], F[1], 'b', linewidth=4)
    plt.plot(P[2], F[2], 'b', linewidth=4)
    plt.plot(P[3], F[3], 'b', linewidth=4)
    plt.plot(P[4], F[4], 'b', linewidth=4)
    plt.plot(P[5], F[5], 'b', linewidth=4)
    plt.plot(P[6], F[6], 'b', linewidth=4)
    plt.plot(P[7], F[7], 'b', linewidth=4)
    plt.plot(P[8], F[8], 'b', linewidth=4)
    plt.plot(P[9], F[9], 'r', linewidth=4)


    #ax.set_xlim(xmin=1, xmax=2E3)
    #ax.set_ylim(ymin=2E-5, ymax=2E1)
    #ax.set_xbound(lower=0.5, upper=2E3)
    #plt.show()
    plt.savefig('electrons_momentum9.png', bbox_inches='tight')