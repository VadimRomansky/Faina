from matplotlib import animation
from pylab import *
import numpy as np

def plot_distributions_4():
    Efile0e = open('../examples_data/gamma0.2_theta0-90/Ee3.dat','r').readline().split()
    Ffile0e= open('../examples_data/gamma0.2_theta0-90/Fs3.dat','r').readline().split()

    Efile0p = open('../examples_data/gamma0.2_theta0-90_protons/Ee3.dat', 'r').readline().split()
    Ffile0p = open('../examples_data/gamma0.2_theta0-90_protons/Fs3.dat', 'r').readline().split()

    me = 9.1E-28
    mp = 1.6E-24
    c = 3E10
    c2 = c*c
    q = 4.84E-10

    Ne = len(Efile0e)
    Np = len(Efile0p)

    Nd = 4

    Ee = np.zeros([Nd,Ne])
    Ep = np.zeros([Nd,Np])
    Fe = np.zeros([Nd,Ne])
    Fp = np.zeros([Nd,Np])

    Ee[0] = open('../examples_data/gamma0.2_theta0-90/Ee3.dat', 'r').readline().split()
    Fe[0] = open('../examples_data/gamma0.2_theta0-90/Fs3.dat', 'r').readline().split()

    Ep[0] = open('../examples_data/gamma0.2_theta0-90_protons/Ee3.dat', 'r').readline().split()
    Fp[0] = open('../examples_data/gamma0.2_theta0-90_protons/Fs3.dat', 'r').readline().split()

    Ee[1] = open('../examples_data/gamma0.3_theta0-90/Ee3.dat', 'r').readline().split()
    Fe[1] = open('../examples_data/gamma0.3_theta0-90/Fs3.dat', 'r').readline().split()

    Ep[1] = open('../examples_data/gamma0.3_theta0-90_protons/Ee3.dat', 'r').readline().split()
    Fp[1] = open('../examples_data/gamma0.3_theta0-90_protons/Fs3.dat', 'r').readline().split()

    Ee[2] = open('../examples_data/gamma0.5_theta0-90/Ee3.dat', 'r').readline().split()
    Fe[2] = open('../examples_data/gamma0.5_theta0-90/Fs3.dat', 'r').readline().split()

    Ep[2] = open('../examples_data/gamma0.5_theta0-90_protons/Ee3.dat', 'r').readline().split()
    Fp[2] = open('../examples_data/gamma0.5_theta0-90_protons/Fs3.dat', 'r').readline().split()

    Ee[3] = open('../examples_data/gamma1.5_theta0-90/Ee3.dat', 'r').readline().split()
    Fe[3] = open('../examples_data/gamma1.5_theta0-90/Fs3.dat', 'r').readline().split()

    Ep[3] = open('../examples_data/gamma1.5_theta0-90_protons/Ee3.dat', 'r').readline().split()
    Fp[3] = open('../examples_data/gamma1.5_theta0-90_protons/Fs3.dat', 'r').readline().split()

    dEe = np.zeros([Nd,Ne])
    dEp = np.zeros([Nd,Np])
    Pe = np.zeros([Nd,Ne])
    Pp = np.zeros([Nd,Np])

    massRatio = np.zeros([2])
    massRatio[0] = me/mp
    massRatio[1] = 1

    totalEenergy = np.zeros([Nd])
    totalPenergy = np.zeros([Nd])
    totalEkineticenergy = np.zeros([Nd])
    totalPkineticenergy = np.zeros([Nd])

    for k in range(Nd):
        for i in range(Ne):
            Ee[k,i] = (float(Ee[k,i])+1)*me*c2
            Fe[k,i] = float(Fe[k,i])
            Pe[k,i] = sqrt(Ee[k,i]*Ee[k,i] - me*me*c2*c2)/c
        for i in range(Ne):
            if i == 0:
                dEe[k,i] = Ee[k,i+1]-Ee[k,i]
            else:
                dEe[k,i] = Ee[k,i] - Ee[k,i-1]

        for i in range(Np):
            Ep[k,i] = (float(Ep[k,i])+1)*mp*c2
            Fp[k,i] = float(Fp[k,i])
            Pp[k,i] = sqrt(Ep[k,i]*Ep[k,i] - mp*mp*c2*c2)/c
        for i in range(Np):
            if i == 0:
                dEp[k,i] = Ep[k,i+1]-Ep[k,i]
            else:
                dEp[k,i] = Ep[k,i] - Ep[k,i-1]

        norme = 0
        for i in range(Ne):
            norme = norme + Fe[k,i]*dEe[k,i]

        normp = 0
        for i in range(Np):
            normp = normp + Fp[k,i] * dEp[k,i]

        Emax = (1.4E8)*me*c2
        for i in range(Ne):
            #F[j][i] = F[j][i]*(E[j,i]+massRatio[j])*(E[j,i]+massRatio[j])/norm[j]
            Fe[k,i] = Fe[k,i]/norme
            totalEenergy+= Fe[k,i]*Ee[k,i]*dEe[k,i]
            totalEkineticenergy += Fe[k,i] * (Ee[k,i] - me*c*c) * dEe[k,i]
            Fe[k,i] = (Fe[k,i]*Pe[k,i]*Pe[k,i]*Pe[k,i]*c2/Ee[k,i])/(me*c)
            #Fe[k,i] = Fe[k,i]*Ee[k,i]*Ee[k,i]
            Pe[k,i] = Pe[k,i]/(me*c)

        for i in range(Np):
            # F[j][i] = F[j][i]*(E[j,i]+massRatio[j])*(E[j,i]+massRatio[j])/norm[j]
            Fp[k,i] = Fp[k,i]/ normp
            totalPenergy += Fp[k,i] * Ep[k,i] * dEp[k,i]
            totalPkineticenergy += Fp[k,i] * (Ep[k,i] - mp*c*c) * dEp[k,i]
            Fp[k,i] = (Fp[k,i]*Pp[k,i]*Pp[k,i]*Pp[k,i]*c2/Ep[k,i])/(mp*c)
            #Fp[k, i] = Fp[k, i] * Ep[k, i] * Ep[k, i]
            Pp[k,i] = Pp[k,i]/(mp*c)

        for i in range(Ne):
            #Fe[k,i] = Fe[k,i]*Ee[k,i]*Ee[k,i]
            Ee[k,i]=Ee[k,i]/(me*c2)

        for i in range(Np):
            #Fp[k,i] = Fp[k,i]*Ep[k,i]*Ep[k,i]
            Ep[k,i] = Ep[k,i]/(mp*c2)

    u = 0.25*0.75*c
    l = 0.1*1.9E17
    B = 0.34

    #F3 = np.zeros([Nd,Ne])
    #for k in range(Nd):
    #    for i in range(Ne):
   #         a = 4*q*q*q*q*B*B*Ee[k,i]*mp*c*c*l/(9*(me**4)*(c**7)*u)

    #        Eeshifted[k,i] = Ee[k,i]/(1 + a)
    #        F3[k,i] = Fe[k,i]*((1+a)*(1+a))*(Eeshifted[k,i]/Ee[k,i])*(Eeshifted[k,i]/Ee[k,i])
    #        #F3[i] = Fe[i]*exp(-Ee[i]/Emax)

    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    #plt.subplots_adjust(bottom=0.12)
    #ax.tick_params(axis='both', which='major', labelsize=10)
    #ax.tick_params(axis='both', which='minor', labelsize=8)
    #ax.set_xlabel('$E_{kin}/m_p c^2$', fontsize=40,fontweight='bold')
    ax.set_xlabel('$p/m_p c$', fontsize=40,fontweight='bold')
    #ax.set_ylabel('$f(E) E^2$', fontsize=40, fontweight='bold')
    ax.set_ylabel('$f_p(p/m_p c) (p/m_p c)^4$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlim([0.05, 20])
    ax.set_ylim([0.001, 1])
    #extraticks=[1,100]
    #plt.xticks(list(plt.xticks()[0]))
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)

    for k in range(Nd):
        #plt.plot(Ep[k], Fp[k], linewidth=4)
        plt.plot(Pp[k], Fp[k], linewidth=4)

    ax.legend([r'v = 0.2c',r'v = 0.3c',r'v = 0.5c',r'v = 0.75c'], fontsize="25")

    plt.savefig('distributions4_protons.png', bbox_inches='tight')

    ax.clear()
    f2 = plt.figure(figsize=[10, 10])
    ax = f2.add_subplot(111)
    # plt.subplots_adjust(bottom=0.12)
    # ax.tick_params(axis='both', which='major', labelsize=10)
    # ax.tick_params(axis='both', which='minor', labelsize=8)
    #ax.set_xlabel('$E_{kin}/m_p c^2$', fontsize=40, fontweight='bold')
    ax.set_xlabel('$p/m_e c$', fontsize=40,fontweight='bold')
    #ax.set_ylabel('$f(E) E^2$', fontsize=40, fontweight='bold')
    ax.set_ylabel('$f_e(p/m_e c) (p/m_e c)^4$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlim([0.1, 1E3])
    ax.set_ylim([1E-4, 10])
    # extraticks=[1,100]
    # plt.xticks(list(plt.xticks()[0]))
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)

    for k in range(Nd):
        #plt.plot(Ee[k], Fe[k], linewidth=4)
        plt.plot(Pe[k], Fe[k], linewidth=4)

    ax.legend([r'v = 0.2c', r'v = 0.3c', r'v = 0.5c', r'v = 0.75c'], fontsize="25")

    plt.savefig('distributions4_electrons.png', bbox_inches='tight')


