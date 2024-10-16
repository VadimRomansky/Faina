from matplotlib import animation
from pylab import *
import numpy as np

def plot_distributions2():
    Efile0e = open('../examples_data/gamma1.5_combined_cutoff/Ee3.dat','r').readline().split()
    Ffile0e= open('../examples_data/gamma1.5_combined_cutoff/Fs3.dat','r').readline().split()

    Efile0p = open('../examples_data/gamma1.5_combined_protons/Ee3.dat', 'r').readline().split()
    Ffile0p = open('../examples_data/gamma1.5_combined_protons/Fs3.dat', 'r').readline().split()

    me = 9.1E-28
    mp = 1.6E-24
    c = 3E10
    c2 = c*c
    q = 4.84E-10

    Ne = len(Efile0e)
    Np = len(Efile0p)

    Ee = np.zeros([Ne])
    Ep = np.zeros([Np])
    Fe = np.zeros([Ne])
    Fp = np.zeros([Np])
    dEe = np.zeros([Ne])
    dEp = np.zeros([Np])
    Pe = np.zeros([Ne])
    Pp = np.zeros([Np])

    massRatio = np.zeros([2])
    massRatio[0] = me/mp
    massRatio[1] = 1

    totalEenergy = 0
    totalPenergy = 0
    totalEkineticenergy = 0
    totalPkineticenergy = 0

    for i in range(Ne):
        Ee[i] = (float(Efile0e[i])+1)*me*c2
        Fe[i] = float(Ffile0e[i])
        Pe[i] = sqrt(Ee[i]*Ee[i] - me*me*c2*c2)/c
    for i in range(Ne):
        if i == 0:
            dEe[i] = Ee[i+1]-Ee[i]
        else:
            dEe[i] = Ee[i] - Ee[i-1]

    for i in range(Np):
        Ep[i] = (float(Efile0p[i])+1)*mp*c2
        Fp[i] = float(Ffile0p[i])
        Pp[i] = sqrt(Ep[i]*Ep[i] - mp*mp*c2*c2)/c
    for i in range(Np):
        if i == 0:
            dEp[i] = Ep[i+1]-Ep[i]
        else:
            dEp[i] = Ep[i] - Ep[i-1]

    norme = 0
    for i in range(Ne):
            norme = norme + Fe[i]*dEe[i]

    normp = 0
    for i in range(Np):
        normp = normp + Fp[i] * dEp[i]

    Emax = (1.4E8)*me*c2
    for i in range(Ne):
            #F[j][i] = F[j][i]*(E[j,i]+massRatio[j])*(E[j,i]+massRatio[j])/norm[j]
            Fe[i] = Fe[i]/norme
            totalEenergy+= Fe[i]*Ee[i]*dEe[i]
            totalEkineticenergy += Fe[i] * (Ee[i] - me*c*c) * dEe[i]
            Fe[i] = (Fe[i]*Pe[i]*Pe[i]*Pe[i]*c2/Ee[i])/(mp*c)
            Pe[i] = Pe[i]/(mp*c)

    for i in range(Np):
        # F[j][i] = F[j][i]*(E[j,i]+massRatio[j])*(E[j,i]+massRatio[j])/norm[j]
        Fp[i] = Fp[i]/ normp
        totalPenergy += Fp[i] * Ep[i] * dEp[i]
        totalPkineticenergy += Fp[i] * (Ep[i] - mp*c*c) * dEp[i]
        Fp[i] = (Fp[i]*Pp[i]*Pp[i]*Pp[i]*c2/Ep[i])/(mp*c)
        Pp[i] = Pp[i]/(mp*c)

    for i in range(Ne):
        Fe[i] = Fe[i]*Ee[i]*Ee[i]
        Ee[i]=Ee[i]/(mp*c2)

    for i in range(Np):
        Fp[i] = Fp[i]*Ep[i]*Ep[i]
        Ep[i] = Ep[i]/(mp*c2)

    Eeshifted = np.zeros([Ne])
    u = 0.25*0.75*c
    l = 0.1*1.9E17
    B = 0.34

    F3 = np.zeros([Ne])
    for i in range(Ne):
        a = 4*q*q*q*q*B*B*Ee[i]*mp*c*c*l/(9*(me**4)*(c**7)*u)

        Eeshifted[i] = Ee[i]/(1 + a)
        F3[i] = Fe[i]*((1+a)*(1+a))*(Eeshifted[i]/Ee[i])*(Eeshifted[i]/Ee[i])
        #F3[i] = Fe[i]*exp(-Ee[i]/Emax)

    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    #plt.subplots_adjust(bottom=0.12)
    #ax.tick_params(axis='both', which='major', labelsize=10)
    #ax.tick_params(axis='both', which='minor', labelsize=8)
    #ax.set_xlabel('$\\rm E_{kin}/m_p c^2$', fontsize=40,fontweight='bold')
    ax.set_xlabel('$\\rm p/m_p c$', fontsize=40,fontweight='bold')
    #ax.set_ylabel('$\\rm f(E) E^2$', fontsize=40, fontweight='bold')
    ax.set_ylabel('$\\rm f(p/m_p c) (p/m_p c)^4$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    ax.set_xscale("log")
    #extraticks=[1,100]
    #plt.xticks(list(plt.xticks()[0]))
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    #ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])
    plt.plot(Eeshifted, F3, 'b', linewidth=4)
    plt.plot(Ep, Fp, 'r', linewidth=4)
    plt.plot(Ee, Fe, 'b', linewidth=4, linestyle='dashed')
    #plt.plot(Pe, F3, 'b', linewidth=4)
    #plt.plot(Pp, Fp, 'r', linewidth=4)
    #plt.plot(Pe, Fe, 'b', linewidth=4, linestyle='dashed')
    P1=Pe[137:150]
    F1=Fe[137:150]
    P2 = Pp[170:190]
    F2 = Fp[170:190]
    #plt.plot(P1, F1, 'b', linewidth=50, linestyle='dashed', dashes=(5, 1), alpha=0.5)
    #plt.plot(E2, F2, 'r', linewidth=50, linestyle='dashed', dashes=(5, 1), alpha=0.5)
    ax.legend([r'electrons',r'protons'], fontsize="25")
    #ax.set_xlim(xmin=0.5E-3, xmax=5E8)
    #ax.set_ylim(ymin=1E-9, ymax=1E-3)
    #ax.set_xticks([1E-2, 1, 1E2, 1E4, 1E6, 1E8])
    #ax.set_xbound(lower=0.5, upper=2E3)
    #plt.show()
    #plt.savefig('distributions2_energy.png', bbox_inches='tight')
    plt.savefig('distributions2_momentum.png', bbox_inches='tight')

    i1 = 0
    i2 = 0

    for i in range(Ne):
        if Ee[i] > 1000*mp*c2:
            i1 = i
            break

    for i in range(Np):
        if Ep[i] > 1000*mp*c2:
            i2 = i
            break

    print('injection ratio = ', Fp[i2]/Fe[i1])

    print('total E energy = ', totalEenergy)
    print('total P energy = ', totalPenergy)
    print('total E kinetic energy = ', totalEkineticenergy)
    print('total P kinetic energy = ', totalPkineticenergy)

    print('electron initial kinetic energy = ', 0.5*me*c*c)
    print('proton initial kinetic energy = ', 0.5*mp*c*c)

    print('epsilon e =', totalEkineticenergy/(0.5*mp*c*c))
    print('epsilon p =', totalPkineticenergy / (0.5 * mp * c * c))