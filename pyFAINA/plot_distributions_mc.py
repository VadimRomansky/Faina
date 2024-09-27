from matplotlib import animation
from pylab import *
import numpy as np

def plot_distributions_mc():
    Efile0 = open('../examples_data/pdf_sf_gamma1.5/GLE_pdf_sf_B0_0003.dat','r')
    Efile1 = open('../examples_data/pdf_sf_gamma1.5/GLE_pdf_sf_B0_00003.dat','r')
    Efile2 = open('../examples_data/pdf_sf_gamma1.5/GLE_pdf_sf_B0_000003.dat','r')

    lines0 = Efile0.readlines()
    lines1 = Efile1.readlines()
    lines2 = Efile2.readlines()
    N = len(lines0)

    E = np.zeros([3,N])
    F = np.zeros([3,N])
    dE = np.zeros([3,N])
    P = np.zeros([3,N])
    dP = np.zeros([3, N])

    mp = 1.67E-24
    me = 0.91E-27
    c = 3*10E10
    c2 = c*c

    m = me

    factor = 1.2

    for i in range(N):
        s = lines0[i].split()

        P[0, i] = 10**float(s[0])
        F[0, i] = float(s[1])/(mp*mp*mp*c*c*c)

        s = lines1[i].split()

        P[1, i] = 10 ** float(s[0])
        F[1, i] = float(s[1])/(mp*mp*mp*c*c*c)

        s = lines2[i].split()

        P[2, i] = 10 ** float(s[0])
        F[2, i] = float(s[1])/(mp*mp*mp*c*c*c)


    for i in range(N):
        for j in range(3):
            if i == 0 :
                dP[j,i] = P[j][1] - P[j][0]
            else:
                dP[j,i] = P[j][i] - P[j][i-1]


    norm = np.zeros([3])
    for i in range(N):
        for j in range(3):
            norm[j] = norm[j] + F[j][i]*dP[j][i]/(P[j][i]*P[j][i])

    for i in range(N):
        for j in range(3):
            F[j][i] = F[j][i]/norm[j]

    #for i in range(N):
    #    for j in range(6):
    #        F[j][i] = F[j][i]*E[j,i]*E[j,i]

    #for i in range(N):
    #    for j in range(6):
    #        F[j][i] = (F[j][i]*P[j][i]*P[j][i]*P[j][i]*c2/E[j][i])/(m*c)
   #         P[j][i] = P[j][i]/(m*c)

    #norm = np.zeros([6])
    #for i in range(N):
    #    for j in range(6):
    #        dP = 0
    #        if i == 0:
    #            dP = P[j][1] - P[j][0]
    #        else:
    #            dP = P[j][i] - P[j][i-1]
    #        ratio1 = dP/dE[j][i]
    #        ratio2 = E[j][i]/(P[j][i]*m*m*c2*c2)
    #        norm[j] = norm[j] + (F[j][i]/(P[j][i]*P[j][i])) * dP

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
    #plt.xticks(list(plt.xticks()[0]))
    #ax.tick_params(axis='x', size=10, width=4)
    #ax.tick_params(axis='y', size=10, width=4)
    #ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])
    plt.plot(P[0], F[0], linewidth=4)
    plt.plot(P[1], F[1], linewidth=4)
    plt.plot(P[2], F[2], linewidth=4)
    ax.legend([r'$B_{up} = 0.0003 G$',
               r'$B_{up} = 0.00003 G$',
               r'$B_{up} = 0.000003 G$'], fontsize="25")
    ax.set_xlim(xmin=0.01, xmax=1E10)
    ax.set_ylim(ymin=1E-3, ymax=2E0)
    #ax.set_xbound(lower=0.5, upper=2E3)
    #plt.show()
    plt.savefig('distribution_mc.png', bbox_inches='tight')