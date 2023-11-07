from matplotlib import animation
from pylab import *
import numpy as np

def plot_distributions2():
    Efile0 = open('../examples_data/gamma0.66_combined/Ee3.dat','r').readline().split()
    Ffile0= open('../examples_data/gamma0.66_combined/Fs3.dat','r').readline().split()

    Efile0a = open('../examples_data/gamma0.66_combined_protons/Ee3.dat', 'r').readline().split()
    Ffile0a = open('../examples_data/gamma0.66_combined_protons/Fs3.dat', 'r').readline().split()

    me = 9.1E-28
    mp = 1.6E-24
    c = 3E10

    Ne = len(Efile0)
    Np = len(Efile0a)

    Ee = np.zeros([Ne])
    Ep = np.zeros([Np])
    Fe = np.zeros([Ne])
    Fp = np.zeros([Np])
    dEe = np.zeros([Ne])
    dEp = np.zeros([Np])

    massRatio = np.zeros([2])
    massRatio[0] = me/mp
    massRatio[1] = 1

    for i in range(Ne):
        Ee[i] = float(Efile0[i])*me/mp
        Fe[i] = float(Ffile0[i])
        if i == 0:
            dEe[i] = Ee[i+1]-Ee[i]
        else:
            dEe[i] = Ee[i] - Ee[i-1]

    for i in range(Np):
        Ep[i] = float(Efile0a[i])
        Fp[i] = float(Ffile0a[i])
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

    Emax = (1.4E8)*me/mp
    for i in range(Ne):
            #F[j][i] = F[j][i]*(E[j,i]+massRatio[j])*(E[j,i]+massRatio[j])/norm[j]
            Fe[i] = Fe[i]*(Ee[i])*(Ee[i])/norme

    for i in range(Np):
        # F[j][i] = F[j][i]*(E[j,i]+massRatio[j])*(E[j,i]+massRatio[j])/norm[j]
        Fp[i] = Fp[i] * (Ep[i]) * (Ep[i]) / normp

    F3 = np.zeros([Ne])
    for i in range(Ne):
        F3[i] = Fe[i]*exp(-Ee[i]/Emax)

    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    #plt.subplots_adjust(bottom=0.12)
    #ax.tick_params(axis='both', which='major', labelsize=10)
    #ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.set_xlabel('$E_{kin}/m_p c^2$', fontsize=40,fontweight='bold')
    ax.set_ylabel('$f(E) E^2$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    ax.set_xscale("log")
    #extraticks=[1,100]
    #plt.xticks(list(plt.xticks()[0]))
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    #ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])
    plt.plot(Ee, F3, 'b', linewidth=4)
    plt.plot(Ep, Fp, 'r', linewidth=4)
    plt.plot(Ee, Fe, 'b', linewidth=4, linestyle='dashed')
    E1=Ee[137:150]
    F1=Fe[137:150]
    E2 = Ep[170:190]
    F2 = Fp[170:190]
    plt.plot(E1, F1, 'b', linewidth=50, linestyle='dashed', dashes=(5, 1), alpha=0.5)
    #plt.plot(E2, F2, 'r', linewidth=50, linestyle='dashed', dashes=(5, 1), alpha=0.5)
    ax.legend([r'electrons',r'protons'], fontsize="25")
    ax.set_xlim(xmin=0.5E-3, xmax=5E8)
    ax.set_ylim(ymin=0.5E-6, ymax=0.5)
    ax.set_xticks([1E-2, 1, 1E2, 1E4, 1E6, 1E8])
    #ax.set_xbound(lower=0.5, upper=2E3)
    #plt.show()
    plt.savefig('distributions.png', bbox_inches='tight')

    i1 = 0
    i2 = 0

    for i in range(Ne):
        if Ee[i] > 1000:
            i1 = i
            break

    for i in range(Np):
        if Ep[i] > 1000:
            i2 = i
            break

    print('injection ratio = ', Fp[i2]/Fe[i1])